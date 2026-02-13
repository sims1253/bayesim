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
#' @param resume Logical; if TRUE, attempt to resume from checkpoint
#' @param force_restart Logical; if TRUE, ignore existing checkpoint and restart
#' @param progress Logical; if TRUE, show progress bar
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
#'   metrics = c("rmse", "bias"),
#'   n_replicates = 100L,
#'   seed = 42L
#' )
#' result <- run_simulation(config)
#' }
run_simulation <- function(
  config,
  resume = FALSE,
  force_restart = FALSE,
  progress = TRUE
) {
  # Validate config
  if (!is_simulation_config(config)) {
    cli::cli_abort("config must be a SimulationConfig object")
  }
  validate_simulation_config(config)

  timer <- make_timer()
  timer$start()

  # Set up global RNG
  setup_global_rng(config@seed)

  # Compute config fingerprint for checkpoint validation
  config_fingerprint <- compute_config_fingerprint(config)

  # Initialize or resume from checkpoint
  should_resume <- resume && !force_restart
  has_existing_checkpoints <- !is.null(config@result_path) &&
    dir.exists(file.path(config@result_path, "checkpoints"))

  if (should_resume && has_existing_checkpoints) {
    resume_data <- tryCatch(
      load_for_resume(config@result_path, config, force_restart = FALSE),
      error = function(e) NULL
    )
    if (!is.null(resume_data)) {
      cli::cli_alert_info("Resuming from previous run")
      task_grid <- resume_data$task_grid
      prior_results <- resume_data$prior_results
      # Don't call init_checkpoint_dir() here!
    } else {
      # Resume failed, start fresh
      cli::cli_alert_info("Resume failed, starting fresh")
      init_checkpoint_dir(config@result_path, config_fingerprint)
      task_grid <- create_task_grid(config)
      prior_results <- NULL
    }
  } else {
    # Fresh start
    if (!is.null(config@result_path)) {
      if (force_restart && can_resume(config@result_path)) {
        cli::cli_alert_info("Force restart: ignoring previous run")
      }
      init_checkpoint_dir(config@result_path, config_fingerprint)
    }
    task_grid <- create_task_grid(config)
    prior_results <- NULL
  }

  total_tasks <- nrow(task_grid)

  # Convert config to spec for worker transport
  config_spec <- as_config_spec(config)

  # Add data_generator to config_spec for run_task
  config_spec$data_generator <- config@data_generator

  # Resolve metrics
  metrics <- resolve_metrics_from_registry(config@metrics)

  cli::cli_alert_info("Starting simulation with {total_tasks} tasks")

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

  # Write final checkpoint with merged results
  if (!is.null(config@result_path)) {
    # Create combined task_results for checkpoint
    write_checkpoint(
      config@result_path,
      results$task_grid,
      results$task_results,
      config_fingerprint
    )
  }

  # Build final result with merged results
  build_simulation_result(
    config = config,
    task_results = results$task_results,
    task_grid = results$task_grid,
    final_results_df = final_results_df,
    elapsed = timer$elapsed(),
    checkpoint_path = config@result_path
  )
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

  task_results <- vector("list", nrow(task_grid))
  error_count <- 0
  checkpoint_counter <- 0

  if (progress && n_pending > 0) {
    pb <- cli::cli_progress_bar(
      total = n_pending,
      format = "Running tasks {cli::pb_current}/{cli::pb_total} [{cli::pb_elapsed}]"
    )
  }

  for (i in seq_len(n_pending)) {
    task <- get_task_spec(
      task_grid,
      pending$task_id[i],
      config
    )

    result <- run_task_safe(
      task = task,
      config_spec = config_spec,
      fitter = fitter,
      metrics = metrics,
      retain = retain
    )

    # Store result with bounds checking
    idx <- match(task$task_id, task_grid$task_id)
    if (is.na(idx)) {
      cli::cli_abort("Task ID '{task$task_id}' not found in task_grid")
    }
    task_results[[idx]] <- result
    task_grid <- update_task_status(task_grid, task$task_id, result$status)

    # Track errors
    if (result$status == "failed") {
      error_count <- error_count + 1
      if (error_count >= max_errors) {
        cli::cli_alert_warning("Reached max_errors ({max_errors}), stopping")
        break
      }
    }

    if (progress) {
      cli::cli_progress_update(id = pb)
    }

    # Periodic checkpointing
    if (!is.null(result_path) && !is.null(config_fingerprint)) {
      checkpoint_counter <- checkpoint_counter + 1
      if (checkpoint_counter >= checkpoint_every) {
        write_checkpoint(
          result_path,
          task_grid,
          task_results,
          config_fingerprint
        )
        checkpoint_counter <- 0
      }
    }
  }

  if (progress) {
    cli::cli_progress_done(id = pb)
  }

  list(
    task_results = task_results,
    task_grid = task_grid
  )
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

  # Extract errors
  errors <- do.call(
    rbind,
    lapply(task_results, function(tr) {
      if (is.null(tr) || is.null(tr$error)) {
        return(NULL)
      }
      data.frame(
        task_id = tr$task_id,
        error_class = tr$error$class,
        error_message = tr$error$message,
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
    summary = summary,
    timing = list(total = elapsed),
    errors = errors,
    checkpoint_path = checkpoint_path
  )
}
