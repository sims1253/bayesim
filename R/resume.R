#' @keywords internal
NULL

#' Check if Resumable Run Exists
#'
#' Checks whether a valid checkpoint exists that can be resumed.
#' A run can be resumed if:
#' - Both run_manifest.json and latest.json exist
#' - latest.json points to a valid checkpoint_id
#' - The referenced checkpoint can be read and validated
#'
#' @param result_path Character; path to results directory containing checkpoints.
#'
#' @return TRUE if a valid run can be resumed, FALSE otherwise.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' if (can_resume("/path/to/results")) {
#'   summary <- get_resume_summary("/path/to/results")
#'   cli::cli_alert_info("Found {summary$n_completed} completed tasks")
#' }
#' }
can_resume <- function(result_path) {
  if (is.null(result_path)) {
    return(FALSE)
  }

  manifest_path <- file.path(result_path, "run_manifest.json")
  latest_path <- file.path(result_path, "latest.json")

  if (!file.exists(manifest_path) || !file.exists(latest_path)) {
    return(FALSE)
  }

  # Check for valid checkpoint reference
  latest <- tryCatch(
    jsonlite::read_json(latest_path),
    error = function(e) NULL
  )

  if (is.null(latest) || is.null(latest$checkpoint_id)) {
    return(FALSE)
  }

  # Verify checkpoint can be read
  checkpoint <- tryCatch(
    read_checkpoint(result_path, latest$checkpoint_id),
    error = function(e) NULL
  )

  !is.null(checkpoint)
}

#' Load Run for Resume
#'
#' Loads previous run state for resumption. Validates schema compatibility
#' and configuration fingerprint before allowing resume.
#'
#' @param result_path Character; path to results directory containing checkpoints.
#' @param config SimulationConfig; current configuration object.
#' @param force_restart Logical; if TRUE, ignore fingerprint mismatch and
#'   restart anyway. Default is FALSE.
#'
#' @return A list with elements:
#'   - `task_grid`: Task grid with restored status from checkpoint
#'   - `prior_results`: Data frame of results from checkpoint
#'   - `checkpoint_id`: ID of the checkpoint being resumed from
#'
#' @details
#' The function performs the following validation steps:
#' 1. Reads and validates run_manifest.json
#' 2. Checks schema version compatibility
#' 3. Computes and compares configuration fingerprint (unless force_restart)
#' 4. Finds the most recent valid checkpoint
#' 5. Rebuilds task grid with status from checkpoint
#'
#' @export
#'
#' @examples
#' \dontrun{
#' config <- simulation_config(...)
#' resume_state <- load_for_resume("/path/to/results", config)
#' task_grid <- resume_state$task_grid
#' prior_results <- resume_state$prior_results
#' }
load_for_resume <- function(result_path, config, force_restart = FALSE) {
  # Read manifest
  manifest_path <- file.path(result_path, "run_manifest.json")
  manifest <- tryCatch(
    jsonlite::read_json(manifest_path),
    error = function(e) {
      cli::cli_abort(
        "Cannot read run manifest: {manifest_path}",
        class = "bayesim_checkpoint_error"
      )
    }
  )

  # Check run schema version compatibility
  if (!identical(manifest$run_schema_version, RUN_SCHEMA_VERSION)) {
    cli::cli_abort(
      c(
        "Incompatible checkpoint schema",
        "x" = "Checkpoint run_schema_version: {manifest$run_schema_version}",
        "i" = "Expected version: {RUN_SCHEMA_VERSION}"
      ),
      class = "bayesim_checkpoint_error"
    )
  }

  # Check result schema version compatibility
  if (!identical(manifest$result_schema_version, RESULT_SCHEMA_VERSION)) {
    cli::cli_abort(
      c(
        "Incompatible result schema",
        "x" = "Checkpoint result_schema_version: {manifest$result_schema_version}",
        "i" = "Expected version: {RESULT_SCHEMA_VERSION}"
      ),
      class = "bayesim_checkpoint_error"
    )
  }

  # Find valid checkpoint (scan backward from latest if corrupted)
  checkpoint <- get_latest_valid_checkpoint(result_path)
  if (is.null(checkpoint)) {
    cli::cli_abort(
      "No valid checkpoint found in {result_path}",
      class = "bayesim_checkpoint_error"
    )
  }

  # Validate checkpoint fingerprint (not just manifest)
  expected_fingerprint <- compute_config_fingerprint(config)
  if (!identical(checkpoint$meta$config_fingerprint, expected_fingerprint)) {
    if (!force_restart) {
      cli::cli_abort(
        c(
          "Configuration fingerprint mismatch",
          "x" = "Cannot resume: checkpoint was created with different configuration",
          "i" = "Use force_restart = TRUE to restart anyway"
        ),
        class = "bayesim_checkpoint_error"
      )
    }
  }

  # Rebuild task grid with restored status
  fresh_grid <- create_task_grid(config)

  # Merge status from checkpoint ledger
  task_grid <- merge_task_grid_status(fresh_grid, checkpoint$task_grid)

  list(
    task_grid = task_grid,
    prior_results = checkpoint$results_df,
    checkpoint_id = checkpoint$checkpoint_id
  )
}

#' Merge Task Grid Status from Checkpoint
#'
#' Updates a freshly-generated task grid with terminal status values
#' from a checkpoint. Pending tasks remain pending; only tasks that
#' were completed (success/failed/skipped) in the checkpoint are updated.
#'
#' @param fresh_grid Task grid tibble from [create_task_grid()].
#' @param checkpoint_grid Task grid tibble from checkpoint.
#'
#' @return Task grid tibble with status merged from checkpoint.
#'
#' @details
#' The function preserves the deterministic task ordering and RNG seeds
#' from the fresh grid, only updating status for tasks that were terminal
#' in the checkpoint. This ensures resumed runs maintain identical
#' reproducibility guarantees.
#'
#' @keywords internal
merge_task_grid_status <- function(fresh_grid, checkpoint_grid) {
  # Start with fresh grid (has all tasks with pending status and rng_seed)
  # Update status for tasks that were terminal in checkpoint

  # Identify terminal tasks in checkpoint
  terminal_statuses <- c("success", "failed", "skipped")
  terminal_mask <- checkpoint_grid$status %in% terminal_statuses

  if (!any(terminal_mask)) {
    # No terminal tasks to merge
    return(fresh_grid)
  }

  terminal_tasks <- checkpoint_grid[terminal_mask, ]

  # Build lookup for fast status retrieval
  status_lookup <- stats::setNames(
    terminal_tasks$status,
    terminal_tasks$task_id
  )

  # Update status for matching task IDs
  task_ids <- fresh_grid$task_id
  for (i in seq_along(task_ids)) {
    task_id <- task_ids[i]
    if (task_id %in% names(status_lookup)) {
      fresh_grid$status[i] <- status_lookup[task_id]
    }
  }

  fresh_grid
}

#' Merge Prior and New Results
#'
#' Combines results from a resumed run with new results from continued
#' execution. Deduplicates by task_id using last-write-wins semantics
#' (new results take precedence over prior results).
#'
#' @param prior_results Data frame of results from checkpoint.
#' @param new_results Data frame of results from new execution.
#'
#' @return Combined data frame with unique task_id values.
#'
#' @details
#' If both inputs are NULL or empty, returns NULL. If only one input
#' has data, returns that input. Otherwise, removes any task_ids from
#' prior_results that appear in new_results, then combines.
#'
#' @keywords internal
merge_results <- function(prior_results, new_results) {
  # Handle empty/NULL cases
  if (is.null(prior_results) || nrow(prior_results) == 0) {
    return(new_results)
  }

  if (is.null(new_results) || nrow(new_results) == 0) {
    return(prior_results)
  }

  # Remove duplicates from prior (keep new - last-write-wins)
  prior_only <- prior_results[
    !prior_results$task_id %in% new_results$task_id,
    ,
    drop = FALSE
  ]

  # Combine: prior-only rows + all new rows
  combined <- rbind(prior_only, new_results)

  # Reset row names for clean output
  rownames(combined) <- NULL

  combined
}

#' Get Resume Summary
#'
#' Returns a summary of what would be resumed from a checkpoint.
#' Useful for informing users about resume state before actually
#' loading and resuming.
#'
#' @param result_path Character; path to results directory containing checkpoints.
#'
#' @return A list with elements:
#'   - `checkpoint_id`: ID of the checkpoint
#'   - `n_total`: Total number of tasks
#'   - `n_completed`: Number of completed (success + failed) tasks
#'   - `n_pending`: Number of pending tasks
#'
#' Returns NULL if no valid resume state exists.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' summary <- get_resume_summary("/path/to/results")
#' if (!is.null(summary)) {
#'   cli::cli_alert_info("Checkpoint: {summary$checkpoint_id}")
#'   cli::cli_alert_info("Completed: {summary$n_completed}/{summary$n_total}")
#'   cli::cli_alert_info("Remaining: {summary$n_pending}")
#' }
#' }
get_resume_summary <- function(result_path) {
  if (!can_resume(result_path)) {
    return(NULL)
  }

  latest_path <- file.path(result_path, "latest.json")
  latest <- tryCatch(
    jsonlite::read_json(latest_path),
    error = function(e) NULL
  )

  if (is.null(latest) || is.null(latest$checkpoint_id)) {
    return(NULL)
  }

  checkpoint <- tryCatch(
    read_checkpoint(result_path, latest$checkpoint_id),
    error = function(e) NULL
  )

  if (is.null(checkpoint)) {
    return(NULL)
  }

  # Extract metadata
  meta <- checkpoint$meta

  if (is.null(meta)) {
    return(NULL)
  }

  list(
    checkpoint_id = checkpoint$checkpoint_id,
    n_total = meta$n_tasks,
    n_completed = meta$n_success + meta$n_failed,
    n_pending = meta$n_pending
  )
}

#' Format Resume Summary for Display
#'
#' Creates a human-readable summary string of the resume state.
#'
#' @param summary A summary list from [get_resume_summary()].
#'
#' @return A character string suitable for display.
#'
#' @keywords internal
format_resume_summary <- function(summary) {
  if (is.null(summary)) {
    return("No resumable state found")
  }

  sprintf(
    "Checkpoint %s: %d/%d tasks completed (%d pending)",
    summary$checkpoint_id,
    summary$n_completed,
    summary$n_total,
    summary$n_pending
  )
}
