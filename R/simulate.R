#' @title Simulation Study Entry Point
#' @description Main functions for running complete simulation studies with
#'   deterministic reproducibility.
#' @name simulate
#' @keywords internal
NULL

#' Run a Simulation Study
#'
#' Executes a complete simulation study with deterministic reproducibility.
#'
#' @param config A SimulationConfig S7 object
#' @param resume Character strategy controlling how an existing `result_path`
#'   is treated: `"auto"` (default) resumes when the path holds a compatible
#'   run and starts fresh only when the path is absent or empty — an existing
#'   run that is incompatible or corrupt aborts rather than being silently
#'   overwritten, as does a non-empty directory without a run manifest;
#'   `"never"` starts fresh and errors if `result_path` already holds a run
#'   or unrelated files; `"must"` resumes and errors when no compatible
#'   checkpoint exists. Only tasks with terminal status
#'   (`"success"`/`"failed"`) are carried over; all other tasks re-run with
#'   their original RNG streams.
#' @param progress Logical; if TRUE, show progress bar
#' @param workers Positive integer, NULL, or "multisession". When non-NULL,
#'   `mirai::daemons(workers)` is set up for the run and torn down on exit —
#'   the simple path for local parallelism. `workers = 1` is genuinely
#'   sequential: no daemons are launched, so package-external S7 fitters and
#'   metrics keep their method dispatch (S7 method tables are registered per
#'   process and do not travel to daemon workers). Parallel execution starts
#'   at `workers >= 2`. Must be NULL when daemons are already set (use
#'   `mirai::daemons()` directly for the advanced/HPC path: remote daemons,
#'   TLS, etc.). Daemons are set before the model bank ships.
#' @param verbose Logical; if TRUE (default), print preflight, lifecycle, and
#'   completion messages. This is independent of `progress`: set
#'   `progress = FALSE` to hide the task bar while keeping run messages, or
#'   `verbose = FALSE` for a quiet programmatic run.
#'
#' @return A bayesim_simulation_result S3 object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' config <- simulation_config(
#'   data_grid = data.frame(n = c(100, 500)),
#'   fit_grid = data.frame(model = "baseline"),
#'   data_generator = my_data_gen,
#'   fitter = my_fitter,
#'   metrics = list(pred_rmse_metric(), pred_bias_metric()),
#'   n_replicates = 100L,
#'   seed = 42L
#' )
#' result <- run_simulation(config)
#' }
run_simulation <- function(
  config,
  resume = c("auto", "never", "must"),
  progress = TRUE,
  workers = NULL,
  verbose = TRUE
) {
  resume <- match.arg(resume)

  # Validate config
  if (!is_simulation_config(config)) {
    cli::cli_abort("config must be a SimulationConfig object")
  }
  validate_simulation_config(config)

  # C2: `workers` convenience argument. When non-NULL, set up mirai daemons for
  # the run and tear them down on exit — but ONLY when no daemons were already
  # set (respect user-managed daemons). Error if both `workers` and existing
  # daemons are present. This happens before the model-bank everywhere() ship.
  # workers = 1 is genuinely sequential: a one-daemon mirai cluster would run
  # tasks in a separate process with no parallelism benefit, and S7 method
  # tables registered by package-external fitters/metrics (S7::method(<-)) do
  # not exist there, so dispatch would break. Leaving daemons unset makes
  # purrr's in_parallel() fall back to ordinary in-process sequential
  # execution. workers >= 2 retains full parallel semantics.
  if (!is.null(workers)) {
    if (isTRUE(mirai::daemons_set())) {
      stop(bayesim_config_error(c(
        "{.arg workers} is non-NULL but mirai daemons are already set.",
        i = "Pass {.code workers = NULL} to use the existing daemons, or call {.code mirai::daemons(0)} first."
      )))
    }
    sequential_workers <- is.numeric(workers) &&
      length(workers) == 1L &&
      isTRUE(workers == 1)
    if (sequential_workers) {
      # Genuinely sequential: no daemons, tasks run in-process.
    } else {
      mirai::daemons(workers)
      on.exit(mirai::daemons(0), add = TRUE)
    }
  }

  timer <- make_timer()
  timer$start()

  # Derive task streams under L'Ecuyer-CMRG without leaking either the RNG kind
  # or .Random.seed into the caller's session after the run returns.
  withr::local_seed(config@seed, .rng_kind = "L'Ecuyer-CMRG")

  # Convert config to spec for hashing/manifest storage/worker transport
  study_spec <- new_study_spec(config)
  run_policy <- new_run_policy(
    config,
    progress = progress,
    workers = workers,
    verbose = verbose
  )
  manifest_spec <- as_config_spec(config)
  config_spec <- manifest_spec
  config_spec$data_generator <- config@data_generator
  config_spec$result_path <- run_policy$result_path
  config_spec$package_name <- utils::packageName()

  # The StudySpec is the only source of scientific identity. Runtime policy
  # changes (path, retention, checkpoint cadence, error budget, adaptive stop)
  # therefore cannot change the fingerprint used for checkpoint compatibility.
  config_fingerprint <- compute_config_fingerprint(study_spec)
  run_store <- new_run_store(
    result_path = run_policy$result_path,
    config_fingerprint = config_fingerprint,
    config_spec = manifest_spec,
    checkpoint_format = run_policy$checkpoint_format,
    keep_checkpoints = run_policy$keep_checkpoints,
    retention_spec = run_policy$retain,
    run_policy_spec = as_run_policy_spec(run_policy)
  )

  # Initialize or resume from checkpoint. Existing non-empty paths are never
  # silently overwritten: a missing/corrupt/incompatible manifest is an
  # ambiguous collision, not a fresh run.
  path_state <- result_path_state(run_policy$result_path)
  if (
    identical(resume, "never") &&
      path_state %in% c("run", "ambiguous", "file")
  ) {
    cli::cli_abort(
      c(
        "Cannot use {.arg resume = 'never'} with an existing result path",
        x = "The path contains a run or unrelated files.",
        i = "Choose a new result_path or use resume = 'auto'/'must'."
      ),
      class = "bayesim_checkpoint_error"
    )
  }

  has_existing_checkpoints <- identical(path_state, "run") &&
    can_resume(run_policy$result_path)
  resume_data <- NULL

  if (identical(resume, "must") && !has_existing_checkpoints) {
    cli::cli_abort("No compatible checkpoint state found to resume")
  }

  if (!identical(resume, "never") && identical(path_state, "run")) {
    resume_attempt <- tryCatch(
      load_for_resume(run_policy$result_path, config, run_store = run_store),
      error = identity
    )

    if (!inherits(resume_attempt, "error")) {
      resume_data <- resume_attempt
    } else {
      stop(resume_attempt)
    }
  }

  if (!is.null(resume_data)) {
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Resuming from previous run")
    }
    task_grid <- resume_data$task_grid
    prior_results <- resume_data$prior_results
    prior_task_results <- resume_data$prior_task_results
  } else {
    if (!is.null(run_policy$result_path)) {
      run_store$initialize()
    }
    task_grid <- create_task_grid(config)
    prior_results <- NULL
    prior_task_results <- NULL
  }

  metrics <- config@metrics %||% list()

  # S6: fresh one-time-warning state for this run, so repeated runs in one
  # session each warn on their own conditions instead of inheriting a stale
  # "already warned" flag.
  .reset_warn_once()

  # F1: condensed prelight one-liner before the run starts. In condensed mode
  # preflight also warns when metrics need capabilities the fitter lacks (R1).
  # If the preflight path itself errors — e.g. a malformed grid list-column —
  # surface it instead of starting the run in silence (R2).
  if (isTRUE(verbose)) {
    tryCatch(
      preflight(config, condensed = TRUE),
      error = function(e) {
        cli::cli_warn(c(
          "Preflight check could not be completed; the run will proceed.",
          i = conditionMessage(e)
        ))
      }
    )
    cli::cli_alert_info(sprintf(
      "Starting simulation with %d tasks",
      nrow(task_grid)
    ))
  }

  # Build the model bank for BrmsFitter (compile once per distinct spec).
  # The bank is stored in a session option (set_model_bank) and pushed to
  # daemons in execute_tasks(); fit() retrieves it via get_model_bank(). This
  # keeps the bank out of the S7 fitter (which would corrupt the config
  # fingerprint).
  if (
    S7::S7_inherits(config@fitter, BrmsFitter) &&
      isTRUE(config@fitter@precompile)
  ) {
    model_bank <- build_model_bank(
      fitter = config@fitter,
      fit_grid = config@fit_grid,
      data_generator = config@data_generator,
      data_spec_template = as.list(config@data_grid[1L, , drop = FALSE]),
      result_path = run_policy$result_path
    )
    set_model_bank(model_bank)
    # F6: clear the session bank after the run so it does not leak across runs
    # (a stale bank from a previous run could mismatch a new fit_grid).
    on.exit(set_model_bank(NULL), add = TRUE)
  } else {
    set_model_bank(NULL)
  }

  # Execute tasks with periodic checkpointing
  results <- execute_tasks(
    task_grid = task_grid,
    config = config,
    config_spec = config_spec,
    fitter = config@fitter,
    metrics = metrics,
    retain = run_policy$retain,
    max_errors = run_policy$max_errors,
    progress = run_policy$progress,
    result_path = run_policy$result_path,
    config_fingerprint = config_fingerprint,
    checkpoint_every = run_policy$checkpoint_every,
    keep_checkpoints = run_policy$keep_checkpoints,
    prior_results_df = prior_results %||%
      data.frame(
        task_id = character(),
        status = character()
      ),
    prior_task_results = prior_task_results,
    adaptive_next_check = resume_data$adaptive_next_check %||% NULL,
    adaptive_state = resume_data$adaptive_state %||% NULL,
    stop_on = run_policy$stop_on,
    verbose = run_policy$verbose,
    run_store = run_store
  )

  timer$stop()

  # Convert new results to dataframe
  new_results_df <- results_to_dataframe(results$task_results)

  # Merge with prior results if resuming
  if (!is.null(prior_results) && nrow(prior_results) > 0) {
    final_results_df <- merge_results(prior_results, new_results_df)
  } else {
    final_results_df <- new_results_df
  }

  final_task_results <- materialize_task_results(
    results$task_results,
    final_results_df,
    results$task_grid,
    prior_task_results = prior_task_results
  )

  # Write final checkpoint with merged results
  if (!is.null(config_fingerprint)) {
    run_store$write(
      task_grid = results$task_grid,
      task_results = final_task_results,
      prior_results_df = data.frame(
        task_id = character(),
        status = character()
      ),
      prior_task_results = final_task_results,
      adaptive_next_check = results$adaptive_next_check %||% NULL,
      adaptive_state = results$adaptive_state %||% NULL
    )
  }

  # Build final result with merged results
  result <- build_simulation_result(
    config = config,
    task_results = final_task_results,
    task_grid = results$task_grid,
    final_results_df = final_results_df,
    elapsed = timer$elapsed(),
    checkpoint_path = run_policy$result_path
  )

  # F7 (issue #53): one end-of-run block — completion status, task counts,
  # stop reason, results path, and the literal resume command when work
  # remains. Subsumes the F2 failure detail. Gated by verbose only.
  if (isTRUE(verbose)) {
    print_run_summary(result)
  }

  result
}

#' Execute all tasks in grid
#'
#' @param task_grid A task grid tibble from create_task_grid()
#' @param config A SimulationConfig S7 object
#' @param config_spec Plain list config spec for worker transport
#' @param fitter S7 Fitter object
#' @param metrics List of Metric objects
#' @param retain Character vector of what to retain
#' @param max_errors Maximum number of errors before stopping
#' @param progress Logical; if TRUE, show progress bar
#' @param result_path Character; path for checkpoint storage (optional)
#' @param config_fingerprint Character; configuration fingerprint for validation
#' @param checkpoint_every Integer; write checkpoint every N completed tasks.
#'   B4: also bounds the number of task results held in memory at once.
#' @param keep_checkpoints Integer; number of checkpoint commit directories
#'   retained. Pruning removes old commit directories only; immutable outcome
#'   shards and ledger history are never pruned.
#' @param prior_results_df Previously resumed result rows, cached in memory to
#'   avoid re-reading and re-hashing the prior checkpoint for every batch.
#' @param run_store Optional internal RunStore adapter used for checkpoint
#'   persistence. Defaults to an adapter created from `result_path`.
#' @param stop_on Optional adaptive stopping policy from the RunPolicy.
#' @param verbose Logical; if TRUE, print lifecycle messages independently of
#'   the task progress bar.
#'
#' @return A list with task_results and task_grid
#'
#' @keywords internal
execute_tasks <- function(
  task_grid,
  config,
  config_spec,
  fitter,
  metrics,
  retain,
  max_errors,
  progress,
  result_path = NULL,
  config_fingerprint = NULL,
  checkpoint_every = 50L,
  keep_checkpoints = 2L,
  prior_results_df = data.frame(task_id = character(), status = character()),
  prior_task_results = NULL,
  adaptive_next_check = NULL,
  adaptive_state = NULL,
  stop_on = NULL,
  verbose = TRUE,
  run_store = NULL
) {
  if (is.null(run_store)) {
    run_store <- new_run_store(
      result_path = result_path,
      config_fingerprint = config_fingerprint,
      checkpoint_format = config@checkpoint_format,
      keep_checkpoints = keep_checkpoints
    )
  }
  if (!"stop_reason" %in% names(task_grid)) {
    task_grid$stop_reason <- NA_character_
  }
  pending <- get_pending_tasks(task_grid)
  n_pending <- nrow(pending)
  pending_indices <- match(pending$task_id, task_grid$task_id)
  if (anyNA(pending_indices)) {
    stop(bayesim_internal_error(
      "Pending task IDs could not be resolved in the task grid"
    ))
  }
  # Adaptive studies proceed in replicate rounds: each condition cell gets
  # replicate 1, then replicate 2, and so on. This makes a precision check
  # representative across cells instead of allowing the first cell to consume
  # all its replicates before other conditions are observed.
  if (!is.null(stop_on) && n_pending > 0L) {
    pending_indices <- pending_indices[order(
      task_grid$rep_idx[pending_indices],
      task_grid$data_idx[pending_indices],
      task_grid$fit_idx[pending_indices]
    )]
  }
  # B4: one knob. batch_size = checkpoint_every (also bounds in-memory results).
  batch_size <- as.integer(checkpoint_every)

  # For memory-bounded execution, we use a list that may contain NULLs
  # for checkpointed results that have been cleared from memory.
  # The final results are assembled from checkpoint data on disk.
  task_results <- vector("list", nrow(task_grid))
  if (!is.null(prior_task_results) && length(prior_task_results) > 0L) {
    prior_outcomes <- Filter(
      function(x) !is.null(x) && is_bayesim_task_result(x),
      prior_task_results
    )
    prior_ids <- vapply(prior_outcomes, function(x) x$task_id, character(1))
    prior_hits <- match(prior_ids, task_grid$task_id)
    keep_prior <- !is.na(prior_hits)
    task_results[prior_hits[keep_prior]] <- prior_outcomes[keep_prior]
  }
  # Failed outcomes are terminal and survive resume, so they have already
  # consumed the run's error budget. Counting only failures from this process
  # would let an unchanged policy execute work that the prior run stopped.
  error_count <- sum(task_grid$status == "failed", na.rm = TRUE)
  error_budget_exhausted <- is.finite(max_errors) &&
    ((max_errors == 0 && error_count > 0) ||
      (max_errors > 0 && error_count >= max_errors))

  # I3: optional adaptive stopping policy (NULL => run all tasks). The check
  # fires after each batch whose completed-task count reaches a check_every
  # boundary AND >= min_reps. bayesim_adaptive_evaluate never throws.
  if (!is.null(stop_on)) {
    adaptive_next_check <- as.integer(
      adaptive_next_check %||% stop_on$min_reps
    )
    if (is.na(adaptive_next_check) || adaptive_next_check < 1L) {
      adaptive_next_check <- stop_on$min_reps
    }
  }
  adaptive_stopped <- FALSE
  tasks_since_checkpoint <- 0L

  # C1: exactly one progress system. purrr's .progress drives the per-task bar
  # inside each batch (passed through run_batch); the outer loop only emits
  # per-batch checkpoint messages.

  if (n_pending > 0 && !error_budget_exhausted) {
    # F6: ship the model bank and run the daemon_setup hook ONCE per
    # execute_tasks() invocation, BEFORE the batch loop. Previously this lived
    # in run_batch() and re-serialized the full bank to all daemons every batch
    # (e.g. 200x for a 10k-task run at batch 50). Note: daemons launched
    # mid-run by the user will not have the bank — daemons must be set before
    # run_simulation().
    if (isTRUE(mirai::daemons_set())) {
      model_bank <- get_model_bank()
      if (!is.null(model_bank)) {
        mirai::everywhere(
          options(bayesim.model_bank = mb),
          .args = list(mb = model_bank)
        )
      }
      daemon_setup <- config@daemon_setup
      if (is.function(daemon_setup)) {
        mirai::everywhere(
          daemon_setup(),
          .args = list(daemon_setup = daemon_setup)
        )
      }
    }

    batch_start <- 1L

    while (batch_start <= n_pending) {
      # A zero error budget means "stop after the first error", not "run one
      # task and stop even when it succeeds". Positive budgets still bound a
      # batch so a checkpoint never overshoots the configured failure budget.
      remaining_error_budget <- if (is.finite(max_errors) && max_errors > 0) {
        max(1L, as.integer(max_errors - error_count))
      } else {
        batch_size
      }

      current_batch_size <- min(batch_size, remaining_error_budget)
      if (is.finite(max_errors) && max_errors == 0) {
        current_batch_size <- 1L
      }
      if (!is.null(stop_on)) {
        completed_before <- completed_replicate_rounds(task_grid)
        rounds_needed <- max(
          1L,
          as.integer(adaptive_next_check - completed_before)
        )
        remaining_indices <- pending_indices[batch_start:n_pending]
        remaining_reps <- task_grid$rep_idx[remaining_indices]
        round_cutoff <- min(remaining_reps) + rounds_needed - 1L
        adaptive_batch_size <- sum(remaining_reps <= round_cutoff)
        current_batch_size <- min(current_batch_size, adaptive_batch_size)
      }
      batch_end <- min(batch_start + current_batch_size - 1L, n_pending)
      batch_indices <- pending_indices[batch_start:batch_end]
      batch_tasks <- lapply(batch_indices, function(row_idx) {
        get_task_spec_at(task_grid, row_idx, config)
      })

      batch_results <- run_batch(
        batch_tasks = batch_tasks,
        config_spec = config_spec,
        fitter = fitter,
        metrics = metrics,
        retain = retain,
        progress = progress
      )

      task_results[batch_indices] <- batch_results
      batch_statuses <- vapply(batch_results, `[[`, character(1), "status")
      task_grid$status[batch_indices] <- batch_statuses
      error_count <- error_count + sum(batch_statuses == "failed")
      tasks_since_checkpoint <- tasks_since_checkpoint + length(batch_indices)

      # C1: re-raise fatal conditions after collecting the batch. run_task_safe
      # captured them (rather than throwing across the daemon boundary) and
      # marked error$fatal with the full condition class chain. Before
      # re-raising, checkpoint the successful and recoverable outcomes this
      # batch produced: a fatal error kills the run, but its siblings are
      # genuine completed outcomes and must survive for resume. The fatal task
      # itself is reset to pending — it never produced a scientific outcome,
      # so resume must retry it rather than treat it as terminal.
      fatal_result <- batch_results[
        vapply(
          batch_results,
          function(r) {
            isTRUE(r$status == "failed") && isTRUE(r$error$fatal)
          },
          logical(1)
        )
      ]
      if (length(fatal_result) > 0L) {
        fr <- fatal_result[[1]]
        fatal_positions <- match(
          vapply(fatal_result, function(r) r$task_id, character(1)),
          task_grid$task_id
        )
        task_grid$status[fatal_positions] <- "pending"
        if ("stop_reason" %in% names(task_grid)) {
          task_grid$stop_reason[fatal_positions] <- NA_character_
        }
        task_results[fatal_positions] <- list(NULL)

        if (!is.null(config_fingerprint)) {
          tryCatch(
            {
              non_null_indices <- which(
                !vapply(task_results, is.null, logical(1))
              )
              run_store$write(
                task_grid = task_grid,
                task_results = task_results[non_null_indices],
                prior_results_df = prior_results_df,
                prior_task_results = prior_task_results,
                adaptive_next_check = adaptive_next_check,
                adaptive_state = adaptive_state
              )
            },
            error = function(e) {
              cli::cli_warn(c(
                "Could not checkpoint completed outcomes before stopping.",
                i = conditionMessage(e)
              ))
              NULL
            }
          )
        }

        err <- fr$error
        cond_class <- err$condition_class %||%
          c("bayesim_internal_error", "bayesim_error", "error", "condition")
        cond_class <- unique(c(cond_class, "error", "condition"))
        stop(structure(
          list(message = err$error_message, call = NULL),
          class = cond_class
        ))
      }

      # I3: adaptive stopping. Check whenever a scheduled threshold has been
      # crossed, rather than relying on exact modulo equality (checkpoint and
      # adaptive intervals need not divide one another). The check runs before
      # lightening so per-task metrics are still in memory.
      if (!is.null(stop_on)) {
        n_completed <- sum(task_grid$status == "success", na.rm = TRUE)
        check_every <- as.integer(stop_on$check_every %||% 50L)
        min_reps <- as.integer(stop_on$min_reps %||% 50L)
        completed_rounds <- completed_replicate_rounds(task_grid)
        while (
          !adaptive_stopped &&
            completed_rounds >= adaptive_next_check
        ) {
          evaluation <- bayesim_adaptive_evaluate(
            task_results,
            task_grid,
            stop_on,
            config,
            completed_rounds = completed_rounds,
            checked_at = n_completed
          )
          stop_now <- evaluation$stop
          adaptive_state <- evaluation$state
          if (isTRUE(stop_now)) {
            if (isTRUE(verbose)) {
              cli::cli_alert_info(
                "Adaptive stop: MCSE of '{stop_on$measure}' for '{stop_on$estimand}' below {stop_on$target_mcse} after {n_completed} completed tasks; skipping remaining tasks"
              )
            }
            # Mark all still-pending tasks with an explicit policy reason. They
            # remain resumable because merge_task_grid_status() only treats
            # success and failed as terminal.
            pending_remaining <- which(task_grid$status == "pending")
            for (ti in pending_remaining) {
              task_grid$status[ti] <- "skipped"
              task_grid$stop_reason[ti] <- "adaptive_stop"
              task_results[[ti]] <- new_task_result(
                task_id = task_grid$task_id[ti],
                status = "skipped",
                timing = list(total = 0),
                stop_reason = "adaptive_stop",
                error = list(
                  error_class = "policy_stop",
                  error_message = "Task not executed (adaptive stop)"
                )
              )
            }
            adaptive_stopped <- TRUE
          }
          adaptive_next_check <- adaptive_next_check + check_every
          adaptive_state$next_check <- adaptive_next_check
        }
      }

      error_budget_exhausted <- is.finite(max_errors) &&
        ((max_errors == 0 && error_count > 0) ||
          (max_errors > 0 && error_count >= max_errors))
      stopping <- adaptive_stopped || error_budget_exhausted
      last_batch <- batch_end >= n_pending
      checkpoint_due <- tasks_since_checkpoint >= batch_size ||
        stopping ||
        last_batch

      if (!is.null(config_fingerprint) && checkpoint_due) {
        non_null_indices <- which(!vapply(task_results, is.null, logical(1)))
        results_to_checkpoint <- task_results[non_null_indices]

        run_store$write(
          task_grid = task_grid,
          task_results = results_to_checkpoint,
          prior_results_df = prior_results_df,
          prior_task_results = prior_task_results,
          adaptive_next_check = adaptive_next_check,
          adaptive_state = adaptive_state
        )
        tasks_since_checkpoint <- 0L

        for (j in non_null_indices) {
          tr <- task_results[[j]]
          if (!is.null(tr)) {
            task_results[[j]] <- lighten_task_result(
              tr,
              retention_for_task_result(retain, tr$status, tr$warnings)
            )
          }
        }
      }

      if (stopping) {
        if (error_budget_exhausted) {
          if (isTRUE(verbose)) {
            cli::cli_alert_warning(
              "Reached max_errors ({max_errors}), stopping"
            )
          }
        }
        break
      }

      batch_start <- batch_end + 1L
    }
  }

  # After the main loop, fill any remaining NULL slots with policy-stopped
  # outcomes. These rows are deliberately not terminal on resume.
  for (i in seq_along(task_results)) {
    if (
      is.null(task_results[[i]]) &&
        is_resumable_task_status(task_grid$status[i])
    ) {
      task_grid$status[i] <- "skipped"
      task_grid$stop_reason[i] <- if (adaptive_stopped) {
        "adaptive_stop"
      } else {
        "max_errors"
      }
      task_results[[i]] <- new_task_result(
        task_id = task_grid$task_id[i],
        status = "skipped",
        timing = list(total = 0),
        stop_reason = task_grid$stop_reason[i],
        error = list(
          error_class = "policy_stop",
          error_message = "Task not executed (run stopped)"
        )
      )
    }
  }

  list(
    task_results = task_results,
    task_grid = task_grid,
    adaptive_next_check = adaptive_next_check,
    adaptive_state = adaptive_state
  )
}

# C1: restore_mirai_condition() and the miraiError/errorValue inspection loop
# were deleted. run_task_safe() is now total (fatal conditions are captured into
# failed task results with the full class chain), so transport carries only
# bayesim_task_result objects; the controller re-raises fatal conditions after
# collecting a batch (see execute_tasks()).

#' Run one batch of tasks
#'
#' Dispatches a batch of simulation tasks via purrr's mirai integration
#' (`purrr::map()` + `purrr::in_parallel()`). With no daemons set, purrr falls
#' back to sequential execution automatically, so there is a single code path
#' (C1). mirai remains the daemon engine; daemons/model bank/`daemon_setup` are
#' managed by [execute_tasks()].
#'
#' `run_task_safe()` is total (never throws), so transport is pure transport:
#' every returned element is a `bayesim_task_result`. Transport-level failures
#' (e.g. daemon death) surface as errors from `purrr::map()` and propagate
#' unchanged out of [execute_tasks()] and `run_simulation()`: they abort the
#' run, and outcomes from the interrupted batch are lost to that call (the last
#' committed checkpoint still holds everything before it).
#'
#' @param batch_tasks List of task specification lists
#' @param config_spec Plain list config spec for worker transport
#' @param fitter S7 Fitter object
#' @param metrics List of Metric objects
#' @param retain Character vector of what to retain
#' @param progress Logical; if TRUE, surface purrr's progress display
#'
#' @return A list of bayesim_task_result objects
#'
#' @keywords internal
run_batch <- function(
  batch_tasks,
  config_spec,
  fitter,
  metrics,
  retain,
  progress = FALSE
) {
  # C1: single dispatch path. in_parallel() crates via carrier::crate(), which
  # strips the lambda's environment; the mapped element is the ONLY positional
  # argument, and every other dependency must be declared as a named constant
  # (purrr in_parallel convention). run_task_safe resolves on daemons because
  # bayesim is installed there (its callee run_task resolves in the bayesim
  # namespace on the daemon). User generators/fitters are crated into the
  # transport and so need not be installed on daemons themselves.
  results <- purrr::map(
    batch_tasks,
    purrr::in_parallel(
      \(task) run_task_safe(task, config_spec, fitter, metrics, retain),
      run_task_safe = run_task_safe,
      config_spec = config_spec,
      fitter = fitter,
      metrics = metrics,
      retain = retain
    ),
    .progress = if (progress) "bayesim tasks" else FALSE
  )
  results
}

materialize_task_results <- function(
  task_results,
  final_results_df,
  task_grid,
  prior_task_results = NULL
) {
  if (is.null(final_results_df) || nrow(final_results_df) == 0) {
    return(task_results)
  }

  result_lookup <- split(final_results_df, final_results_df$task_id)
  prior_lookup <- if (!is.null(prior_task_results)) {
    prior_task_results <- Filter(
      function(x) !is.null(x) && is_bayesim_task_result(x),
      prior_task_results
    )
    stats::setNames(
      prior_task_results,
      vapply(prior_task_results, function(x) x$task_id, character(1))
    )
  } else {
    list()
  }

  lapply(seq_along(task_results), function(i) {
    if (!is.null(task_results[[i]])) {
      return(task_results[[i]])
    }

    canonical <- prior_lookup[[task_grid$task_id[i]]]
    if (!is.null(canonical)) {
      return(canonical)
    }

    row <- result_lookup[[task_grid$task_id[i]]]
    if (is.null(row) || nrow(row) == 0) {
      return(NULL)
    }

    # Legacy fallback (flat summary rows only): rebuild the canonical outcome
    # via the shared migration helper so truth, stop_reason, error, and
    # diagnostics are routed to their proper fields instead of $metrics.
    task_outcomes_from_dataframe(row[1, , drop = FALSE])[[1]]
  })
}

#' Build final simulation result
#'
#' @param config A SimulationConfig S7 object
#' @param task_results List of bayesim_task_result objects
#' @param task_grid The task grid tibble with updated statuses
#' @param final_results_df Dataframe with merged results (prior + new)
#' @param elapsed Total elapsed time in seconds
#' @param checkpoint_path Character; path where checkpoints were stored
#'
#' @return A bayesim_simulation_result S3 object
#'
#' @keywords internal
build_simulation_result <- function(
  config,
  task_results,
  task_grid,
  final_results_df = NULL,
  elapsed,
  checkpoint_path = NULL
) {
  # Use merged results if available, otherwise compute from task_results
  if (!is.null(final_results_df) && nrow(final_results_df) > 0) {
    summary <- final_results_df
  } else {
    # Flatten task results to summary tibble
    summary_rows <- lapply(task_results, function(tr) {
      if (is.null(tr)) {
        return(NULL)
      }
      row <- list(
        task_id = tr$task_id,
        status = tr$status,
        total_time = tr$timing$total
      )
      if (!is.null(tr$metrics)) {
        row <- c(row, tr$metrics)
      }
      if (!is.null(tr$diagnostics)) {
        row <- c(row, tr$diagnostics)
      }
      row
    })

    summary <- tibble::as_tibble(
      do.call(rbind, lapply(summary_rows, as.data.frame))
    )
  }

  # Enrich summary with data_grid, fit_grid, and rep_idx columns
  summary <- enrich_summary_with_grid_columns(
    summary = summary,
    task_grid = task_grid,
    data_grid = config@data_grid,
    fit_grid = config@fit_grid
  )

  # Extract errors
  errors <- do.call(
    rbind,
    lapply(task_results, function(tr) {
      if (is.null(tr) || is.null(tr$error)) {
        return(NULL)
      }
      data.frame(
        task_id = tr$task_id,
        error_class = tr$error$error_class,
        error_message = tr$error$error_message,
        stringsAsFactors = FALSE
      )
    })
  )

  if (is.null(errors)) {
    errors <- data.frame(
      task_id = character(),
      error_class = character(),
      error_message = character(),
      stringsAsFactors = FALSE
    )
  }

  # I8: optional parquet sidecar for the summary. The rds checkpoint remains
  # the canonical resume artifact; this parquet file is for downstream
  # consumption (pandas/arrow/polars). Best-effort: warn on failure.
  if (
    identical(config@summary_format, "parquet") &&
      !is.null(checkpoint_path)
  ) {
    summary_parquet_path <- file.path(checkpoint_path, "summary.parquet")
    tryCatch(
      write_results_parquet(summary, summary_parquet_path),
      error = function(e) {
        cli::cli_warn(c(
          "Failed to write parquet summary to {.file {summary_parquet_path}}.",
          i = conditionMessage(e)
        ))
      }
    )
  }

  # E4: record each metric's declared summary_type so summarize_simulation()
  # can pick the right aggregation/MCSE without name heuristics.
  metric_summary_types <- NULL
  metric_field_metadata <- NULL
  metrics <- config@metrics %||% list()
  if (length(metrics) > 0) {
    metric_summary_types <- lapply(metrics, function(m) {
      if (S7::S7_inherits(m)) m@summary_type else "mean"
    })
    names(metric_summary_types) <- vapply(
      metrics,
      function(m) {
        if (S7::S7_inherits(m)) m@name else "unknown"
      },
      character(1)
    )
    metric_field_metadata <- lapply(metrics, function(m) {
      if (S7::S7_inherits(m, Metric)) {
        m@schema
      } else {
        list()
      }
    })
    names(metric_field_metadata) <- names(metric_summary_types)
  }

  new_simulation_result(
    config_fingerprint = compute_config_fingerprint(new_study_spec(config)),
    task_results = task_results,
    task_grid = task_grid,
    summary = summary,
    timing = list(total = elapsed),
    errors = errors,
    checkpoint_path = checkpoint_path,
    metric_summary_types = metric_summary_types,
    metric_field_metadata = metric_field_metadata
  )
}

#' Enrich summary with grid columns
#'
#' Adds data_grid columns (prefixed with "data_"), fit_grid columns (prefixed
#' with "fit_"), and rep_idx to the summary tibble.
#'
#' @param summary A tibble with task results
#' @param task_grid The task grid tibble with data_idx, fit_idx, rep_idx columns
#' @param data_grid A data frame with data configuration rows
#' @param fit_grid A data frame with fit configuration rows
#'
#' @return The summary tibble with additional columns
#'
#' @keywords internal
enrich_summary_with_grid_columns <- function(
  summary,
  task_grid,
  data_grid,
  fit_grid
) {
  # Handle empty summary
  if (is.null(summary) || nrow(summary) == 0) {
    return(summary)
  }

  # Match summary rows to task_grid by task_id
  # Use match to get indices for each task
  task_ids <- summary$task_id
  grid_indices <- match(task_ids, task_grid$task_id)

  # Get data_idx and fit_idx for each task
  data_indices <- task_grid$data_idx[grid_indices]
  fit_indices <- task_grid$fit_idx[grid_indices]
  rep_indices <- task_grid$rep_idx[grid_indices]

  # Add rep_idx column
  summary$rep_idx <- rep_indices

  # Add data_grid columns with "data_" prefix
  if (!is.null(data_grid) && nrow(data_grid) > 0) {
    data_colnames <- names(data_grid)
    for (col_name in data_colnames) {
      new_col_name <- paste0("data_", col_name)
      summary[[new_col_name]] <- data_grid[[col_name]][data_indices]
    }
  }

  if (is.null(data_grid) && "data_spec" %in% names(task_grid)) {
    data_specs <- task_grid$data_spec[grid_indices]
    data_names <- unique(unlist(lapply(data_specs, names), use.names = FALSE))
    for (col_name in data_names) {
      new_col_name <- paste0("data_", col_name)
      # H: preserve atomic types (numeric stays numeric); NA only for non-scalar.
      vals <- lapply(data_specs, function(spec) {
        value <- spec[[col_name]]
        if (is.null(value) || length(value) != 1 || is.list(value)) {
          NA
        } else {
          value
        }
      })
      summary[[new_col_name]] <- simplify_to_atomic(vals)
    }
  }

  # Add fit_grid columns with "fit_" prefix
  if (!is.null(fit_grid) && nrow(fit_grid) > 0) {
    fit_colnames <- names(fit_grid)
    for (col_name in fit_colnames) {
      new_col_name <- paste0("fit_", col_name)
      summary[[new_col_name]] <- fit_grid[[col_name]][fit_indices]
    }
  }

  if (is.null(fit_grid) && "fit_spec" %in% names(task_grid)) {
    fit_specs <- task_grid$fit_spec[grid_indices]
    fit_names <- unique(unlist(lapply(fit_specs, names), use.names = FALSE))
    for (col_name in fit_names) {
      new_col_name <- paste0("fit_", col_name)
      vals <- lapply(fit_specs, function(spec) {
        value <- spec[[col_name]]
        if (is.null(value) || length(value) != 1 || is.list(value)) {
          NA
        } else {
          value
        }
      })
      summary[[new_col_name]] <- simplify_to_atomic(vals)
    }
  }

  summary
}

# H: coerce a list of scalar values to the simplest atomic vector that holds
# them, preserving numeric/integer/logical types; fall back to character only
# when the values are genuinely non-numeric. NA-safe.
simplify_to_atomic <- function(vals) {
  if (length(vals) == 0L) {
    return(character(0))
  }
  non_na <- vals[!is.na(vals)]
  if (length(non_na) == 0L || all(vapply(non_na, is.numeric, logical(1)))) {
    unlist(vals)
  } else if (all(vapply(non_na, is.logical, logical(1)))) {
    unlist(vals)
  } else {
    vapply(
      vals,
      function(v) if (is.na(v)) NA_character_ else as.character(v),
      character(1)
    )
  }
}

#' Resume a simulation from an existing result directory
#'
#' Calls [run_simulation()] with `resume = "must"`.
#'
#' @param result_path Character path to an existing result directory
#' @param config Optional SimulationConfig object. When omitted, bayesim
#'   rebuilds the config from the run manifest. That works only when every
#'   executable component (data generator, fitter, metrics) is a namespaced
#'   package function or class. The run manifest cannot rehydrate
#'   script-defined closures, such as the return value of
#'   [fixed_truth_generator()] or any generator defined in your script: R can
#'   serialize closures, but configless resume intentionally does not restore
#'   arbitrary executable closures. Those runs require the original `config`.
#' @param progress Logical; if TRUE, show progress bar
#' @param verbose Logical; if TRUE, print lifecycle messages independently of
#'   the progress bar.
#'
#' @return A bayesim_simulation_result S3 object
#' @export
resume_simulation <- function(
  result_path,
  config = NULL,
  progress = TRUE,
  verbose = TRUE
) {
  if (is.null(config)) {
    config <- rehydrate_config_from_manifest(result_path)
  }

  if (!is_simulation_config(config)) {
    cli::cli_abort("config must be a SimulationConfig object")
  }

  config@result_path <- result_path
  run_simulation(
    config,
    resume = "must",
    progress = progress,
    verbose = verbose
  )
}

# Workstream I3: adaptive stopping -----------------------------------------

#' Build a quick summary tibble from in-memory task results (I3)
#'
#' Internal helper for adaptive stopping: flattens the non-NULL entries of
#' `task_results` to a summary data.frame via [results_to_dataframe()] and
#' enriches it with data_grid/fit_grid/rep_idx columns (matching what
#' [build_simulation_result()] produces). Used by the internal adaptive
#' evaluator so it can call [performance_measures()] mid-run.
#'
#' @param task_results List of `bayesim_task_result` (possibly with NULLs).
#' @param task_grid The task grid tibble (with up-to-date statuses).
#' @param config A SimulationConfig.
#' @return A data.frame summary, or NULL if no non-NULL results.
#' @keywords internal
bayesim_adaptive_summary <- function(task_results, task_grid, config) {
  non_null <- which(!vapply(task_results, is.null, logical(1)))
  if (length(non_null) == 0L) {
    return(NULL)
  }
  df <- results_to_dataframe(task_results[non_null])
  if (is.null(df) || nrow(df) == 0L) {
    return(NULL)
  }
  df <- enrich_summary_with_grid_columns(
    summary = df,
    task_grid = task_grid,
    data_grid = config@data_grid,
    fit_grid = config@fit_grid
  )
  grid_indices <- match(df$task_id, task_grid$task_id)
  df$.data_idx <- task_grid$data_idx[grid_indices]
  df$.fit_idx <- task_grid$fit_idx[grid_indices]
  df
}

#' Minimum successful replicate round completed in every condition cell.
#'
#' @param task_grid Current task grid with up-to-date statuses.
#' @return A non-negative integer replicate-round count.
#' @keywords internal
completed_replicate_rounds <- function(task_grid) {
  success <- task_grid[task_grid$status == "success", , drop = FALSE]
  if (nrow(success) == 0L) {
    return(0L)
  }
  cells <- interaction(
    task_grid$data_idx,
    task_grid$fit_idx,
    drop = TRUE,
    lex.order = TRUE
  )
  success_cells <- interaction(
    success$data_idx,
    success$fit_idx,
    drop = TRUE,
    lex.order = TRUE
  )
  counts <- table(factor(success_cells, levels = levels(cells)))
  as.integer(min(counts))
}

#' Evaluate adaptive precision and return a persistable decision snapshot.
#' @keywords internal
bayesim_adaptive_evaluate <- function(
  task_results_so_far,
  task_grid,
  stop_on,
  config,
  completed_rounds = completed_replicate_rounds(task_grid),
  checked_at = sum(task_grid$status == "success", na.rm = TRUE)
) {
  if (is.null(stop_on)) {
    return(list(stop = FALSE, state = NULL))
  }
  tryCatch(
    {
      df <- bayesim_adaptive_summary(task_results_so_far, task_grid, config)
      if (is.null(df) || nrow(df) == 0L) {
        return(list(stop = FALSE, state = NULL))
      }
      pm <- performance_measures(
        df,
        estimand = stop_on$estimand,
        by = c(".data_idx", ".fit_idx")
      )
      if (is.null(pm) || nrow(pm) == 0L) {
        return(list(stop = FALSE, state = NULL))
      }
      sel <- pm$measure == stop_on$measure
      if (!any(sel)) {
        return(list(stop = FALSE, state = NULL))
      }
      selected <- pm[sel, , drop = FALSE]
      expected_cells <- nrow(unique(task_grid[c("data_idx", "fit_idx")]))
      enough <- !(nrow(selected) != expected_cells ||
        any(selected$n_sim < stop_on$min_reps))
      mcse <- selected$mcse
      target <- stop_on$target_mcse
      stop <- isTRUE(enough) &&
        length(mcse) > 0L &&
        all(is.finite(mcse) & mcse < target)
      cells <- lapply(seq_len(nrow(selected)), function(i) {
        list(
          data_idx = if (".data_idx" %in% names(selected)) {
            selected$.data_idx[[i]]
          } else {
            i
          },
          fit_idx = if (".fit_idx" %in% names(selected)) {
            selected$.fit_idx[[i]]
          } else {
            1L
          },
          n_sim = selected$n_sim[[i]],
          mcse = selected$mcse[[i]]
        )
      })
      list(
        stop = stop,
        state = list(
          triggered = stop,
          estimand = stop_on$estimand,
          measure = stop_on$measure,
          target_mcse = stop_on$target_mcse,
          checked_at = as.integer(checked_at),
          completed_rounds = as.integer(completed_rounds),
          cells = cells
        )
      )
    },
    error = function(e) {
      # Keep the non-throwing contract (a broken precision check must not kill
      # the run), but do not swallow the failure silently: surface it once per
      # run with the condition message. The detail is also persisted in
      # state$cells$error for post-hoc inspection.
      .warn_once(
        "adaptive_evaluate_error",
        c(
          "Adaptive precision evaluation failed; continuing without a stop decision.",
          i = conditionMessage(e)
        )
      )
      list(
        stop = FALSE,
        state = list(
          triggered = FALSE,
          estimand = stop_on$estimand,
          measure = stop_on$measure,
          target_mcse = stop_on$target_mcse,
          checked_at = as.integer(checked_at),
          completed_rounds = as.integer(completed_rounds),
          cells = list(),
          error = conditionMessage(e)
        )
      )
    }
  )
}

#' Should the run stop early based on the adaptive-stopping policy?
#'
#' @inheritParams bayesim_adaptive_evaluate
#' @return `TRUE` when every condition cell meets the precision target.
#' @keywords internal
bayesim_adaptive_check <- function(
  task_results_so_far,
  task_grid,
  stop_on,
  config
) {
  isTRUE(
    bayesim_adaptive_evaluate(
      task_results_so_far,
      task_grid,
      stop_on,
      config
    )$stop
  )
}
