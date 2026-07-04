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
#' @param resume Character strategy: "auto", "never", or "must"
#' @param progress Logical; if TRUE, show progress bar
#' @param workers Positive integer, NULL, or "multisession". When non-NULL,
#'   `mirai::daemons(workers)` is set up for the run and torn down on exit —
#'   the simple path for local parallelism. Must be NULL when daemons are
#'   already set (use `mirai::daemons()` directly for the advanced/HPC path:
#'   remote daemons, TLS, etc.). Daemons are set before the model bank ships.
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
#'   metrics = list(rmse_metric(), bias_metric()),
#'   n_replicates = 100L,
#'   seed = 42L
#' )
#' result <- run_simulation(config)
#' }
run_simulation <- function(
  config,
  resume = c("auto", "never", "must"),
  progress = TRUE,
  workers = NULL
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
  if (!is.null(workers)) {
    if (isTRUE(mirai::daemons_set())) {
      stop(bayesim_config_error(c(
        "{.arg workers} is non-NULL but mirai daemons are already set.",
        i = "Pass {.code workers = NULL} to use the existing daemons, or call {.code mirai::daemons(0)} first."
      )))
    }
    mirai::daemons(workers)
    on.exit(mirai::daemons(0), add = TRUE)
  }

  timer <- make_timer()
  timer$start()

  # Set up global RNG
  setup_global_rng(config@seed)

  # Convert config to spec for hashing/manifest storage/worker transport
  manifest_spec <- as_config_spec(config)
  config_spec <- manifest_spec
  config_spec$data_generator <- config@data_generator
  config_spec$result_path <- config@result_path
  config_spec$package_name <- utils::packageName()

  # Compute config fingerprint for checkpoint validation
  config_fingerprint <- compute_config_fingerprint(config)

  # Initialize or resume from checkpoint
  has_existing_checkpoints <- !is.null(config@result_path) &&
    can_resume(config@result_path)
  resume_data <- NULL

  if (identical(resume, "must") && !has_existing_checkpoints) {
    cli::cli_abort("No compatible checkpoint state found to resume")
  }

  if (!identical(resume, "never") && has_existing_checkpoints) {
    resume_attempt <- tryCatch(
      load_for_resume(config@result_path, config),
      error = identity
    )

    if (!inherits(resume_attempt, "error")) {
      resume_data <- resume_attempt
    } else if (identical(resume, "must")) {
      stop(resume_attempt)
    } else {
      cli::cli_alert_info(
        "Existing checkpoint state not resumable; starting fresh"
      )
    }
  }

  if (!is.null(resume_data)) {
    cli::cli_alert_info("Resuming from previous run")
    task_grid <- resume_data$task_grid
    prior_results <- resume_data$prior_results
  } else {
    if (!is.null(config@result_path)) {
      init_checkpoint_dir(
        config@result_path,
        config_fingerprint,
        config_spec = manifest_spec,
        checkpoint_format = config@checkpoint_format
      )
    }
    task_grid <- create_task_grid(config)
    prior_results <- NULL
  }

  metrics <- config@metrics %||% list()

  # F1: condensed prelight one-liner before the run starts.
  tryCatch(preflight(config, condensed = TRUE), error = function(e) NULL)

  cli::cli_alert_info(sprintf(
    "Starting simulation with %d tasks",
    nrow(task_grid)
  ))

  # Build the model bank for BrmsFitter (compile once per distinct spec).
  # The bank is stored in a session option (set_model_bank) and pushed to
  # daemons in execute_tasks(); fit() retrieves it via get_model_bank(). This
  # keeps the bank out of the S7 fitter (which would corrupt the config
  # fingerprint).
  if (inherits(config@fitter, "BrmsFitter") && isTRUE(config@fitter@precompile)) {
    model_bank <- build_model_bank(
      fitter = config@fitter,
      fit_grid = config@fit_grid,
      data_generator = config@data_generator,
      data_spec_template = as.list(config@data_grid[1L, , drop = FALSE]),
      result_path = config@result_path,
      seed = config@seed
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
    retain = config@retain,
    max_errors = config@max_errors,
    progress = progress,
    result_path = config@result_path,
    config_fingerprint = config_fingerprint,
    checkpoint_every = config@checkpoint_every
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
    results$task_grid
  )

  # Write final checkpoint with merged results
  if (!is.null(config@result_path)) {
    write_checkpoint(
      config@result_path,
      results$task_grid,
      final_task_results,
      config_fingerprint,
      checkpoint_format = config@checkpoint_format
    )
  }

  # Build final result with merged results
  result <- build_simulation_result(
    config = config,
    task_results = final_task_results,
    task_grid = results$task_grid,
    final_results_df = final_results_df,
    elapsed = timer$elapsed(),
    checkpoint_path = config@result_path
  )

  # F2: compact failure summary when any task failed.
  if (!is.null(result$summary) && any(result$summary$status == "failed", na.rm = TRUE)) {
    print_failure_summary(result)
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
#' @param checkpoint_every Integer; write checkpoint every N completed tasks
#' @param chunk_size Integer; maximum task results to keep in memory before
#'   forcing a checkpoint write and clearing memory (default: same as checkpoint_every)
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
  checkpoint_every = 50L
) {
  pending <- get_pending_tasks(task_grid)
  n_pending <- nrow(pending)
  # B4: one knob. batch_size = checkpoint_every (also bounds in-memory results).
  batch_size <- as.integer(checkpoint_every)

  # For memory-bounded execution, we use a list that may contain NULLs
  # for checkpointed results that have been cleared from memory.
  # The final results are assembled from checkpoint data on disk.
  task_results <- vector("list", nrow(task_grid))
  error_count <- 0
  n_in_memory <- 0

  # C1: exactly one progress system. purrr's .progress drives the per-task bar
  # inside each batch (passed through run_batch); the outer loop only emits
  # per-batch checkpoint messages.

  if (n_pending > 0) {
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
        mirai::everywhere(daemon_setup(), .args = list(daemon_setup = daemon_setup))
      }
    }

    batch_start <- 1L

    while (batch_start <= n_pending) {
      remaining_error_budget <- if (is.finite(max_errors)) {
        max(1L, as.integer(max_errors - error_count))
      } else {
        batch_size
      }

      current_batch_size <- min(batch_size, remaining_error_budget)
      batch_end <- min(batch_start + current_batch_size - 1L, n_pending)
      batch_ids <- pending$task_id[batch_start:batch_end]
      batch_tasks <- lapply(batch_ids, function(task_id) {
        get_task_spec(task_grid, task_id, config)
      })

      batch_results <- run_batch(
        batch_tasks = batch_tasks,
        config_spec = config_spec,
        fitter = fitter,
        metrics = metrics,
        retain = retain,
        progress = progress
      )

      for (k in seq_along(batch_results)) {
        task <- batch_tasks[[k]]
        result <- batch_results[[k]]
        idx <- match(task$task_id, task_grid$task_id)

        if (is.na(idx)) {
          cli::cli_abort("Task ID '{task$task_id}' not found in task_grid")
        }

        task_results[[idx]] <- result
        task_grid <- update_task_status(task_grid, task$task_id, result$status)
        n_in_memory <- n_in_memory + 1L

        if (identical(result$status, "failed")) {
          error_count <- error_count + 1L
        }
      }

      # C1: re-raise fatal conditions after collecting the batch. run_task_safe
      # captured them (rather than throwing across the daemon boundary) and
      # marked error$fatal with the full condition class chain.
      fatal_result <- batch_results[
        vapply(batch_results, function(r) {
          isTRUE(r$status == "failed") && isTRUE(r$error$fatal)
        }, logical(1))
      ]
      if (length(fatal_result) > 0L) {
        fr <- fatal_result[[1]]
        err <- fr$error
        cond_class <- err$condition_class %||%
          c("bayesim_internal_error", "bayesim_error", "error", "condition")
        cond_class <- unique(c(cond_class, "error", "condition"))
        stop(structure(
          list(message = err$error_message, call = NULL),
          class = cond_class
        ))
      }

      if (!is.null(result_path) && !is.null(config_fingerprint)) {
        non_null_indices <- which(!vapply(task_results, is.null, logical(1)))
        results_to_checkpoint <- task_results[non_null_indices]

        write_checkpoint(
          result_path,
          task_grid,
          results_to_checkpoint,
          config_fingerprint,
          checkpoint_format = config@checkpoint_format
        )

        for (j in non_null_indices) {
          tr <- task_results[[j]]
          if (!is.null(tr)) {
            task_results[[j]] <- lighten_task_result(
              tr,
              retention_for_task_result(retain, tr$status, tr$warnings)
            )
          }
        }

        n_in_memory <- sum(!vapply(task_results, is.null, logical(1)))
        gc(verbose = FALSE)
      }

      if (error_count >= max_errors) {
        cli::cli_alert_warning("Reached max_errors ({max_errors}), stopping")
        break
      }

      batch_start <- batch_end + 1L
    }
  }

  # After the main loop, fill any remaining NULL slots with skipped results
  for (i in seq_along(task_results)) {
    if (
      is.null(task_results[[i]]) && identical(task_grid$status[i], "pending")
    ) {
      task_results[[i]] <- new_task_result(
        task_id = task_grid$task_id[i],
        status = "skipped",
        timing = list(total = 0),
        error = list(
          error_class = "skipped",
          error_message = "Task not executed (max_errors reached)"
        )
      )
      task_grid$status[i] <- "skipped"
    }
  }

  list(
    task_results = task_results,
    task_grid = task_grid
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
#' (e.g. daemon death) surface as errors from `purrr::map()` and are re-raised
#' as a bayesim fatal error by the caller.
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
  task_grid
) {
  if (is.null(final_results_df) || nrow(final_results_df) == 0) {
    return(task_results)
  }

  result_lookup <- split(final_results_df, final_results_df$task_id)

  lapply(seq_along(task_results), function(i) {
    if (!is.null(task_results[[i]])) {
      return(task_results[[i]])
    }

    row <- result_lookup[[task_grid$task_id[i]]]
    if (is.null(row) || nrow(row) == 0) {
      return(NULL)
    }

    row <- row[1, , drop = FALSE]
    excluded <- c(
      "task_id",
      "status",
      "error_class",
      "error_message",
      "timing_total"
    )
    metric_cols <- setdiff(names(row), excluded)
    metrics <- lapply(metric_cols, function(col) row[[col]][[1]])
    names(metrics) <- metric_cols

    if (identical(row$status[[1]], "failed")) {
      new_task_result(
        task_id = row$task_id[[1]],
        status = "failed",
        metrics = NULL,
        timing = list(total = row$timing_total[[1]] %||% 0),
        error = list(
          error_class = row$error_class[[1]] %||% "unknown",
          error_message = row$error_message[[1]] %||% "Task failed"
        )
      )
    } else {
      new_task_result(
        task_id = row$task_id[[1]],
        status = row$status[[1]],
        metrics = metrics,
        timing = list(total = row$timing_total[[1]] %||% 0),
        warnings = character()
      )
    }
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

  new_simulation_result(
    config_fingerprint = compute_config_fingerprint(config),
    task_results = task_results,
    task_grid = task_grid,
    summary = summary,
    timing = list(total = elapsed),
    errors = errors,
    checkpoint_path = checkpoint_path
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
      summary[[new_col_name]] <- vapply(
        data_specs,
        function(spec) {
          value <- spec[[col_name]]
          if (is.null(value) || length(value) != 1 || is.list(value)) {
            NA
          } else {
            as.character(value)
          }
        },
        character(1)
      )
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
      summary[[new_col_name]] <- vapply(
        fit_specs,
        function(spec) {
          value <- spec[[col_name]]
          if (is.null(value) || length(value) != 1 || is.list(value)) {
            NA
          } else {
            as.character(value)
          }
        },
        character(1)
      )
    }
  }

  summary
}

#' Resume a simulation from an existing result directory
#'
#' @param result_path Character path to an existing result directory
#' @param config Optional SimulationConfig object
#' @param progress Logical; if TRUE, show progress bar
#'
#' @return A bayesim_simulation_result S3 object
#' @export
resume_simulation <- function(result_path, config = NULL, progress = TRUE) {
  if (is.null(config)) {
    config <- rehydrate_config_from_manifest(result_path)
  }

  if (!is_simulation_config(config)) {
    cli::cli_abort("config must be a SimulationConfig object")
  }

  config@result_path <- result_path
  run_simulation(config, resume = "must", progress = progress)
}
