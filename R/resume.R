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
#' @keywords internal
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
#' 3. Computes and compares configuration fingerprints
#' 4. Finds the most recent valid checkpoint
#' 5. Rebuilds task grid with status from checkpoint
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' config <- simulation_config(...)
#' resume_state <- load_for_resume("/path/to/results", config)
#' task_grid <- resume_state$task_grid
#' prior_results <- resume_state$prior_results
#' }
load_for_resume <- function(result_path, config) {
  # Read manifest
  manifest_path <- file.path(result_path, "run_manifest.json")
  manifest <- tryCatch(
    jsonlite::read_json(manifest_path),
    error = function(e) {
      stop(bayesim_checkpoint_error(
        paste0("Cannot read run manifest: ", manifest_path)
      ))
    }
  )

  # Check run schema version compatibility
  if (!identical(manifest$run_schema_version, RUN_SCHEMA_VERSION)) {
    stop(bayesim_checkpoint_error(paste0(
      "Incompatible checkpoint schema: run_schema_version ",
      manifest$run_schema_version, ", expected ", RUN_SCHEMA_VERSION
    )))
  }

  # Check result schema version compatibility
  if (!identical(manifest$result_schema_version, RESULT_SCHEMA_VERSION)) {
    stop(bayesim_checkpoint_error(paste0(
      "Incompatible result schema: result_schema_version ",
      manifest$result_schema_version, ", expected ", RESULT_SCHEMA_VERSION
    )))
  }

  manifest_format <- manifest$checkpoint_format %||% "rds"
  if (!identical(manifest_format, config@checkpoint_format)) {
    stop(bayesim_checkpoint_error(paste0(
      "Checkpoint format mismatch: checkpoint uses '", manifest_format,
      "' but config requests '", config@checkpoint_format, "'"
    )))
  }

  # Find valid checkpoint (scan backward from latest if corrupted)
  expected_fingerprint <- compute_config_fingerprint(config)
  checkpoint <- get_latest_valid_checkpoint(
    result_path,
    config_fingerprint = expected_fingerprint
  )
  if (is.null(checkpoint)) {
    stop(bayesim_checkpoint_error(
      paste0("No valid compatible checkpoint found in ", result_path)
    ))
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
#' execution. Duplicate task rows must be identical or an integrity error
#' is raised.
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

  duplicate_ids <- intersect(prior_results$task_id, new_results$task_id)

  if (length(duplicate_ids) > 0) {
    for (task_id in duplicate_ids) {
      prior_row <- prior_results[
        prior_results$task_id == task_id,
        ,
        drop = FALSE
      ]
      new_row <- new_results[new_results$task_id == task_id, , drop = FALSE]

      if (nrow(prior_row) != 1 || nrow(new_row) != 1) {
        stop(bayesim_checkpoint_error(
          paste0("Duplicate terminal rows detected for task_id '", task_id, "'")
        ))
      }

      if (
        !identical(
          normalize_result_row(prior_row),
          normalize_result_row(new_row)
        )
      ) {
        stop(bayesim_checkpoint_error(
          paste(
            "Conflicting duplicate terminal rows detected for task_id",
            shQuote(task_id)
          )
        ))
      }
    }
  }

  prior_only <- prior_results[
    !prior_results$task_id %in% duplicate_ids,
    ,
    drop = FALSE
  ]

  # Combine: prior-only rows + all new rows
  combined <- rbind(prior_only, new_results)

  # Reset row names for clean output
  rownames(combined) <- NULL

  combined
}

normalize_result_row <- function(x) {
  cols <- sort(names(x))
  x <- x[, cols, drop = FALSE]

  lapply(x, function(col) {
    if (is.factor(col)) {
      as.character(col)
    } else {
      col
    }
  })
}

rehydrate_config_from_manifest <- function(result_path) {
  manifest <- read_run_manifest(result_path)

  if (is.null(manifest) || is.null(manifest$config_spec)) {
    stop(bayesim_checkpoint_error(
      "Run manifest does not contain a rehydratable configuration; please supply config explicitly"
    ))
  }

  spec <- manifest$config_spec
  data_generator <- rehydrate_function_spec(spec$data_generator_spec)
  fitter <- rehydrate_s7_spec(spec$fitter_spec)
  metrics <- lapply(spec$metrics_spec %||% list(), rehydrate_s7_spec)

  data_grid <- normalize_manifest_df(spec$data_grid)
  fit_grid <- normalize_manifest_df(spec$fit_grid)
  task_grid <- normalize_manifest_df(spec$task_grid)
  retain <- normalize_manifest_retain(spec$retain)
  max_errors <- normalize_manifest_numeric(spec$max_errors, default = Inf)
  seed <- as.integer(normalize_manifest_numeric(spec$seed, default = 1L))
  n_replicates <- as.integer(normalize_manifest_numeric(
    spec$n_replicates,
    default = 1L
  ))

  simulation_config(
    data_grid = data_grid,
    fit_grid = fit_grid,
    task_grid = task_grid,
    data_generator = data_generator,
    fitter = fitter,
    metrics = metrics,
    n_replicates = n_replicates,
    seed = seed,
    result_path = result_path,
    checkpoint_format = spec$checkpoint_format %||% "rds",
    retain = retain,
    max_errors = max_errors
  )
}

rehydrate_function_spec <- function(spec) {
  if (is.null(spec) || !isTRUE(spec$rehydratable) || is.null(spec$reference)) {
    stop(bayesim_checkpoint_error(
      "Run manifest references non-rehydratable executable components; please supply config explicitly"
    ))
  }

  validate_namespace_version(spec$reference$package, spec$reference$version)

  get(spec$reference$name, envir = asNamespace(spec$reference$package))
}

rehydrate_s7_spec <- function(spec) {
  if (is.null(spec) || (length(spec) == 1 && is.na(spec))) {
    return(NULL)
  }

  if (!isTRUE(spec$rehydratable) || is.null(spec$class)) {
    stop(bayesim_checkpoint_error(
      "Run manifest references non-rehydratable executable components; please supply config explicitly"
    ))
  }

  class_parts <- strsplit(spec$class, "::", fixed = TRUE)[[1]]
  if (length(class_parts) != 2) {
    stop(bayesim_checkpoint_error(
      "Cannot rehydrate S7 object without namespaced class reference"
    ))
  }

  validate_namespace_version(class_parts[[1]], spec$package_version)

  constructor <- get(class_parts[[2]], envir = asNamespace(class_parts[[1]]))
  props <- spec$properties %||% list()
  do.call(constructor, props)
}

validate_namespace_version <- function(package, expected_version = NULL) {
  if (is.null(expected_version)) {
    return(invisible(TRUE))
  }

  current_version <- as.character(getNamespaceVersion(package))
  if (!identical(current_version, expected_version)) {
    stop(bayesim_checkpoint_error(paste0(
      "Component version mismatch during resume rehydration: ",
      "Package '", package, "' version is '", current_version,
      "', expected '", expected_version, "'"
    )))
  }

  invisible(TRUE)
}

normalize_manifest_df <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.list(x) && length(x) == 0) {
    return(NULL)
  }

  if (is.data.frame(x)) {
    return(x)
  }

  if (is.list(x) && all(vapply(x, is.list, logical(1)))) {
    return(bind_rows_safe(x))
  }

  tryCatch(as.data.frame(x), error = function(e) x)
}

normalize_manifest_retain <- function(x) {
  if (is.null(x) || (is.list(x) && length(x) == 0)) {
    return(c("metrics", "diagnostics"))
  }

  if (is.character(x)) {
    return(x)
  }

  if (!is.list(x)) {
    return(c("metrics", "diagnostics"))
  }

  lapply(x, function(value) {
    if (is.character(value)) {
      value
    } else if (is.list(value)) {
      unlist(value, use.names = FALSE)
    } else {
      as.character(value)
    }
  })
}

normalize_manifest_numeric <- function(x, default = NULL) {
  if (is.null(x) || (is.list(x) && length(x) == 0)) {
    return(default)
  }

  if (is.character(x) && identical(x, "Inf")) {
    return(Inf)
  }

  as.numeric(x)
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
#' @keywords internal
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
