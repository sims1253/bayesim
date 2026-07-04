#' @title S3 Result Constructors
#' @description Constructors, validators, and helpers for bayesim result objects.
#'   These S3 classes provide consistent interfaces for handling results from
#'   Bayesian model fitting, task execution, and simulation runs.
#' @name results
#' @keywords internal
NULL

# =============================================================================
# bayesim_fit_result
# =============================================================================

#' Check if object is a bayesim_fit_result
#'
#' @param x Object to check
#' @return `TRUE` if `x` inherits from `"bayesim_fit_result"`, `FALSE` otherwise
#' @keywords internal
is_bayesim_fit_result <- function(x) {
  inherits(x, "bayesim_fit_result")
}

#' Validate a bayesim_fit_result object
#'
#' Performs consistency checks on a bayesim_fit_result object.
#'
#' @param x A bayesim_fit_result object to validate
#' @return The input object, invisibly, if validation passes
#'
#' @section Errors:
#' Throws an error if validation fails, with a message indicating the specific
#' validation problem (e.g., class mismatch, success/error inconsistency,
#' invalid timing, missing draws colnames, etc.).
#'
#' @keywords internal
validate_bayesim_fit_result <- function(x) {
  if (!is_bayesim_fit_result(x)) {
    stop(bayesim_contract_error("Object must have class 'bayesim_fit_result'"))
  }

  if (isTRUE(x$success)) {
    if (!is.null(x$error)) {
      stop(bayesim_contract_error("When success is TRUE, error must be NULL"))
    }
  } else {
    if (is.null(x$error)) {
      stop(bayesim_contract_error(
        "When success is FALSE, error must be non-NULL"
      ))
    }
  }

  if (!is.list(x$timing)) {
    stop(bayesim_contract_error("timing must be a list"))
  }
  if (!is.numeric(x$timing$total) || length(x$timing$total) != 1) {
    stop(bayesim_contract_error("timing$total must be a scalar numeric"))
  }
  if (x$timing$total < 0) {
    stop(bayesim_contract_error("timing$total must be >= 0"))
  }

  if (!is.null(x$draws)) {
    if (!is.matrix(x$draws)) {
      stop(bayesim_contract_error("draws must be a matrix when not NULL"))
    }
    if (is.null(colnames(x$draws))) {
      stop(bayesim_contract_error("draws matrix must have colnames"))
    }
  }

  if (!is.character(x$warnings)) {
    stop(bayesim_contract_error("warnings must be a character vector"))
  }

  if (!is.list(x$diagnostics)) {
    stop(bayesim_contract_error("diagnostics must be a list"))
  }

  invisible(x)
}

#' Create a new bayesim_fit_result object
#'
#' Constructs a result object from a Bayesian model fitting operation.
#'
#' @param success Logical scalar indicating if the fit succeeded
#' @param fit The backend fit object (may be NULL if removed by retention policy)
#' @param draws Matrix of posterior draws (S x P), with column names for parameters
#' @param diagnostics Named list of diagnostic values (scalars or named numeric vectors)
#' @param timing List containing total, warmup, and sample timing in seconds
#' @param warnings Character vector of warning messages captured during fitting
#' @param error A condition object if the fit failed, NULL otherwise
#' @param data_bundle Optional list containing the training/test data used for fitting
#'
#' @return A validated `bayesim_fit_result` object
#'
#' @details
#' The bayesim_fit_result class encapsulates all outputs from a Bayesian model
#' fitting operation, including the posterior draws, diagnostics, timing information,
#' and any warnings or errors encountered.
#'
#' Validation rules:
#' - If `success` is FALSE, `error` must be non-NULL
#' - If `success` is TRUE, `error` must be NULL
#' - `timing$total` must be non-negative
#' - If `draws` is not NULL, it must be a matrix with column names
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' # Successful fit
#' draws <- matrix(rnorm(1000), ncol = 2, nrow = 500)
#' colnames(draws) <- c("alpha", "beta")
#' result <- new_fit_result(
#'   success = TRUE,
#'   draws = draws,
#'   diagnostics = list(rhat = c(alpha = 1.01, beta = 1.00)),
#'   timing = list(total = 10.5, warmup = 5.0, sample = 5.5)
#' )
#'
#' # Failed fit
#' result <- new_fit_result(
#'   success = FALSE,
#'   error = simpleError("Convergence failed"),
#'   diagnostics = list(),
#'   timing = list(total = 2.0, warmup = 2.0, sample = 0)
#' )
#' }
new_fit_result <- function(
  success = TRUE,
  fit = NULL,
  draws = NULL,
  diagnostics = list(),
  timing = list(total = 0, warmup = 0, sample = 0),
  warnings = character(),
  error = NULL,
  data_bundle = NULL
) {
  # Ensure warnings is character vector
  if (is.null(warnings)) {
    warnings <- character()
  }

  # Ensure diagnostics is a list
  if (is.null(diagnostics)) {
    diagnostics <- list()
  }

  # Ensure timing components exist
  if (is.null(timing)) {
    timing <- list(total = 0, warmup = 0, sample = 0)
  }
  if (is.null(timing$total)) {
    timing$total <- 0
  }
  if (is.null(timing$warmup)) {
    timing$warmup <- 0
  }
  if (is.null(timing$sample)) {
    timing$sample <- 0
  }

  result <- structure(
    list(
      success = success,
      fit = fit,
      draws = draws,
      diagnostics = diagnostics,
      timing = timing,
      warnings = warnings,
      error = error,
      data_bundle = data_bundle
    ),
    class = "bayesim_fit_result"
  )

  validate_bayesim_fit_result(result)
}

# =============================================================================
# bayesim_task_result
# =============================================================================

#' Check if object is a bayesim_task_result
#'
#' @param x Object to check
#' @return `TRUE` if `x` inherits from `"bayesim_task_result"`, `FALSE` otherwise
#' @keywords internal
is_bayesim_task_result <- function(x) {
  inherits(x, "bayesim_task_result")
}

#' Validate a bayesim_task_result object
#'
#' Performs consistency checks on a bayesim_task_result object.
#'
#' @param x A bayesim_task_result object to validate
#' @return The input object, invisibly, if validation passes
#'
#' @section Errors:
#' Throws an error if validation fails, with a message indicating the specific
#' validation problem (e.g., class mismatch, invalid task_id or status,
#' missing metrics for successful tasks, missing error for failed tasks, etc.).
#'
#' @keywords internal
validate_bayesim_task_result <- function(x) {
  if (!is_bayesim_task_result(x)) {
    stop(bayesim_contract_error("Object must have class 'bayesim_task_result'"))
  }

  if (!is.character(x$task_id) || length(x$task_id) != 1) {
    stop(bayesim_contract_error("task_id must be a scalar character"))
  }

  valid_statuses <- c("success", "failed", "skipped")
  if (!is.character(x$status) || length(x$status) != 1) {
    stop(bayesim_contract_error("status must be a scalar character"))
  }
  if (!(x$status %in% valid_statuses)) {
    stop(bayesim_contract_error(
      "status must be one of: ",
      paste(valid_statuses, collapse = ", ")
    ))
  }

  if (x$status == "success") {
    if (is.null(x$metrics)) {
      stop(bayesim_contract_error(
        "When status is 'success', metrics must not be NULL"
      ))
    }
  }
  if (x$status == "failed") {
    if (is.null(x$error)) {
      stop(bayesim_contract_error(
        "When status is 'failed', error must not be NULL"
      ))
    }
  }

  if (!is.list(x$timing)) {
    stop(bayesim_contract_error("timing must be a list"))
  }
  if (!is.numeric(x$timing$total) || length(x$timing$total) != 1) {
    stop(bayesim_contract_error("timing$total must be a scalar numeric"))
  }
  if (x$timing$total < 0) {
    stop(bayesim_contract_error("timing$total must be >= 0"))
  }

  if (!is.character(x$warnings)) {
    stop(bayesim_contract_error("warnings must be a character vector"))
  }

  invisible(x)
}

#' Create a new bayesim_task_result object
#'
#' Constructs a result object from a single simulation task execution.
#'
#' @param task_id Character scalar identifying the task
#' @param status Character scalar: one of "success", "failed", or "skipped"
#' @param metrics Named list of computed metrics (NULL if task failed or skipped)
#' @param diagnostics Named list of diagnostic values (NULL if task failed)
#' @param timing List containing timing information, must include `total`
#' @param error NULL, or a list with `error_class` and `error_message` if failed
#' @param warnings Character vector of warning messages
#'
#' @return A validated `bayesim_task_result` object
#'
#' @details
#' The bayesim_task_result class captures the outcome of a single simulation task,
#' including its computed metrics, any diagnostics, timing information, and
#' errors or warnings.
#'
#' Validation rules:
#' - `status` must be one of "success", "failed", or "skipped"
#' - If `status` is "success", `metrics` must not be NULL
#' - If `status` is "failed", `error` must not be NULL
#' - `timing$total` must be non-negative
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' # Successful task
#' result <- new_task_result(
#'   task_id = "task_001",
#'   status = "success",
#'   metrics = list(rmse = 0.05, bias = 0.01),
#'   diagnostics = list(n_eff = 500, rhat = 1.01),
#'   timing = list(total = 5.2)
#' )
#'
#' # Failed task
#' result <- new_task_result(
#'   task_id = "task_002",
#'   status = "failed",
#'   error = list(error_class = "convergence_error", error_message = "R-hat > 1.1"),
#'   timing = list(total = 2.0)
#' )
#' }
new_task_result <- function(
  task_id = character(1),
  status = "success",
  metrics = NULL,
  diagnostics = NULL,
  timing = list(total = 0),
  error = NULL,
  warnings = character(),
  truth = NULL
) {
  # Ensure warnings is character vector
  if (is.null(warnings)) {
    warnings <- character()
  }

  # Ensure timing is a list with total
  if (is.null(timing)) {
    timing <- list(total = 0)
  }
  if (is.null(timing$total)) {
    timing$total <- 0
  }

  result <- structure(
    list(
      task_id = task_id,
      status = status,
      metrics = metrics,
      diagnostics = diagnostics,
      timing = timing,
      error = error,
      warnings = warnings,
      # E1: data-generating truth (named numeric or NULL), always retained and
      # flattened to truth__<param> summary columns for recovery analysis.
      truth = truth
    ),
    class = "bayesim_task_result"
  )

  validate_bayesim_task_result(result)
}

# =============================================================================
# bayesim_simulation_result
# =============================================================================

#' Check if object is a bayesim_simulation_result
#'
#' @param x Object to check
#' @return `TRUE` if `x` inherits from `"bayesim_simulation_result"`, `FALSE` otherwise
#' @keywords internal
is_bayesim_simulation_result <- function(x) {
  inherits(x, "bayesim_simulation_result")
}

#' Validate a bayesim_simulation_result object
#'
#' Performs consistency checks on a bayesim_simulation_result object.
#'
#' @param x A bayesim_simulation_result object to validate
#' @return The input object, invisibly, if validation passes
#'
#' @section Errors:
#' Throws an error if validation fails, with a message indicating the specific
#' validation problem (e.g., class mismatch, invalid config_fingerprint,
#' invalid task_results elements, non-data.frame summary/errors, etc.).
#'
#' @keywords internal
validate_bayesim_simulation_result <- function(x) {
  if (!is_bayesim_simulation_result(x)) {
    stop(bayesim_contract_error(
      "Object must have class 'bayesim_simulation_result'"
    ))
  }

  if (
    !is.character(x$config_fingerprint) || length(x$config_fingerprint) != 1
  ) {
    stop(bayesim_contract_error(
      "config_fingerprint must be a scalar character"
    ))
  }

  if (!is.list(x$task_results)) {
    stop(bayesim_contract_error("task_results must be a list"))
  }
  for (i in seq_along(x$task_results)) {
    if (!is_bayesim_task_result(x$task_results[[i]])) {
      stop(bayesim_contract_error(
        "All elements of task_results must be bayesim_task_result objects"
      ))
    }
  }

  if (!is.null(x$task_grid)) {
    if (!is.data.frame(x$task_grid)) {
      stop(bayesim_contract_error(
        "task_grid must be NULL or a data.frame/tibble"
      ))
    }
  }

  if (!is.data.frame(x$summary)) {
    stop(bayesim_contract_error("summary must be a data.frame or tibble"))
  }

  if (!is.list(x$timing)) {
    stop(bayesim_contract_error("timing must be a list"))
  }
  if (!is.numeric(x$timing$total) || length(x$timing$total) != 1) {
    stop(bayesim_contract_error("timing$total must be a scalar numeric"))
  }
  if (x$timing$total < 0) {
    stop(bayesim_contract_error("timing$total must be >= 0"))
  }

  if (!is.data.frame(x$errors)) {
    stop(bayesim_contract_error("errors must be a data.frame or tibble"))
  }

  if (!is.null(x$checkpoint_path)) {
    if (
      !is.character(x$checkpoint_path) ||
        length(x$checkpoint_path) != 1
    ) {
      stop(bayesim_contract_error(
        "checkpoint_path must be NULL or a scalar character"
      ))
    }
  }

  invisible(x)
}

#' Create a new bayesim_simulation_result object
#'
#' Constructs a result object from a complete simulation run.
#'
#' @param config_fingerprint Character hash uniquely identifying the configuration
#' @param task_results List of `bayesim_task_result` objects
#' @param task_grid Tibble with task grid information (task_id, data_idx, fit_idx, rep_idx, status)
#' @param summary Tibble with one row per task, columns for metrics and diagnostics
#' @param timing List containing total, by_phase, and other timing breakdowns
#' @param errors Tibble of failed tasks with error details
#' @param checkpoint_path Path to the final checkpoint file, or NULL
#'
#' @return A validated `bayesim_simulation_result` object
#'
#' @details
#' The bayesim_simulation_result class encapsulates the complete results of a
#' simulation study, including all individual task results, an aggregated summary,
#' timing breakdowns, error information, and checkpoint location.
#'
#' Validation rules:
#' - `config_fingerprint` must be a scalar character
#' - `task_results` must be a list where all elements are `bayesim_task_result`
#' - `summary` and `errors` must be data frames
#' - `timing$total` must be non-negative
#' - `checkpoint_path` must be NULL or a scalar character
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' # Create a simulation result
#' task1 <- new_task_result(
#'   task_id = "task_001",
#'   status = "success",
#'   metrics = list(rmse = 0.05),
#'   timing = list(total = 5.0)
#' )
#'
#' result <- new_simulation_result(
#'   config_fingerprint = "abc123",
#'   task_results = list(task1),
#'   task_grid = tibble::tibble(
#'     task_id = "task_001",
#'     data_idx = 1L,
#'     fit_idx = 1L,
#'     rep_idx = 1L,
#'     status = "success"
#'   ),
#'   summary = tibble::tibble(task_id = "task_001", rmse = 0.05),
#'   timing = list(total = 10.0, by_phase = list(setup = 2.0, fit = 8.0)),
#'   errors = tibble::tibble(task_id = character(), error_message = character()),
#'   checkpoint_path = "/path/to/checkpoint.rds"
#' )
#' }
new_simulation_result <- function(
  config_fingerprint = character(1),
  task_results = list(),
  task_grid = NULL,
  summary = NULL,
  timing = list(total = 0),
  errors = NULL,
  checkpoint_path = NULL
) {
  # Ensure task_results is a list
  if (is.null(task_results)) {
    task_results <- list()
  }

  # Ensure summary is a data.frame (create empty if NULL)
  if (is.null(summary)) {
    summary <- data.frame()
  }

  # Ensure timing is a list with total
  if (is.null(timing)) {
    timing <- list(total = 0)
  }
  if (is.null(timing$total)) {
    timing$total <- 0
  }

  # Ensure errors is a data.frame (create empty if NULL)
  if (is.null(errors)) {
    errors <- data.frame()
  }

  result <- structure(
    list(
      config_fingerprint = config_fingerprint,
      task_results = task_results,
      task_grid = task_grid,
      summary = summary,
      timing = timing,
      errors = errors,
      checkpoint_path = checkpoint_path
    ),
    class = "bayesim_simulation_result"
  )

  validate_bayesim_simulation_result(result)
}

# =============================================================================
# Print methods
# =============================================================================

#' @exportS3Method
#' @keywords internal
print.bayesim_fit_result <- function(x, ...) {
  cat("<bayesim_fit_result>\n")
  cat("  Success:", if (isTRUE(x$success)) "TRUE" else "FALSE", "\n")
  if (!is.null(x$draws)) {
    cat(
      "  Draws:",
      nrow(x$draws),
      "iterations x",
      ncol(x$draws),
      "parameters\n"
    )
  } else {
    cat("  Draws: NULL\n")
  }
  if (length(x$diagnostics) > 0) {
    cat("  Diagnostics:", length(x$diagnostics), "metrics\n")
  }
  cat("  Total time:", round(x$timing$total, 2), "s\n")
  if (length(x$warnings) > 0) {
    cat("  Warnings:", length(x$warnings), "\n")
  }
  if (!is.null(x$error)) {
    cat("  Error:", conditionMessage(x$error), "\n")
  }
  invisible(x)
}

#' @exportS3Method
#' @keywords internal
print.bayesim_task_result <- function(x, ...) {
  cat("<bayesim_task_result>\n")
  cat("  Task ID:", x$task_id, "\n")
  cat("  Status:", x$status, "\n")
  if (!is.null(x$metrics)) {
    cat("  Metrics:", paste(names(x$metrics), collapse = ", "), "\n")
  }
  if (!is.null(x$diagnostics)) {
    cat("  Diagnostics:", paste(names(x$diagnostics), collapse = ", "), "\n")
  }
  cat("  Total time:", round(x$timing$total, 2), "s\n")
  if (length(x$warnings) > 0) {
    cat("  Warnings:", length(x$warnings), "\n")
  }
  if (!is.null(x$error)) {
    cat("  Error:", x$error$error_message, "\n")
  }
  invisible(x)
}

#' @exportS3Method
#' @keywords internal
print.bayesim_simulation_result <- function(x, ...) {
  cat("<bayesim_simulation_result>\n")
  cat("  Config fingerprint:", x$config_fingerprint, "\n")
  cat("  Tasks:", length(x$task_results), "\n")
  n_success <- sum(vapply(
    x$task_results,
    function(t) t$status == "success",
    logical(1)
  ))
  n_failed <- sum(vapply(
    x$task_results,
    function(t) t$status == "failed",
    logical(1)
  ))
  n_skipped <- sum(vapply(
    x$task_results,
    function(t) t$status == "skipped",
    logical(1)
  ))
  cat("    - Success:", n_success, "\n")
  cat("    - Failed:", n_failed, "\n")
  cat("    - Skipped:", n_skipped, "\n")
  if (!is.null(x$task_grid) && nrow(x$task_grid) > 0) {
    cat(
      "  Task grid:",
      nrow(x$task_grid),
      "rows x",
      ncol(x$task_grid),
      "cols\n"
    )
  }
  cat("  Total time:", round(x$timing$total, 2), "s\n")
  if (!is.null(x$checkpoint_path)) {
    cat("  Checkpoint:", x$checkpoint_path, "\n")
  }
  invisible(x)
}
