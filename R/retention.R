#' @title Memory Management and Retention Policies
#' @description Functions for controlling which artifacts are retained in simulation
#'   results to manage memory usage during large-scale simulations.
#' @name retention
#' @keywords internal
"_PACKAGE"

## Constants ------------------------------------------------------------------

#' Valid retention options
#'
#' Character vector of all valid retention option names that can be used
#' in retention specifications.
#'
#' @format Character vector
#' @keywords internal
RETAIN_OPTIONS <- c(
  "metrics", # Always retained
  "diagnostics", # Convergence diagnostics
  "draws", # Posterior draws matrix
  "predictions", # Predicted values
  "fit", # Raw fit object
  "data", # Input data
  "warnings" # Warning messages
)

#' Retention profiles
#'
#' Pre-defined retention profiles for common use cases. Each profile is a
#' character vector of retention options.
#'
#' @format Named list of character vectors
#' @keywords internal
RETENTION_PROFILES <- list(
  minimal = c("metrics"),
  standard = c("metrics", "diagnostics"),
  debug = c(
    "metrics",
    "diagnostics",
    "draws",
    "predictions",
    "fit",
    "data",
    "warnings"
  )
)

## Resolution ------------------------------------------------------------------

#' Resolve retention specification
#'
#' Converts retention profile name or explicit vector to canonical form.
#' If a profile name is provided, returns the corresponding options.
#' If explicit options are provided, validates them and returns only valid ones.
#'
#' @param retain Character vector of retention options or a profile name
#'   (one of "minimal", "standard", or "debug")
#'
#' @return Character vector of valid retention options
#'
#' @export
#'
#' @examples
#' resolve_retention("minimal")
#' resolve_retention(c("metrics", "draws"))
resolve_retention <- function(retain) {
  # Check if retain is a profile name
  if (
    is.character(retain) &&
      length(retain) == 1 &&
      retain %in% names(RETENTION_PROFILES)
  ) {
    return(RETENTION_PROFILES[[retain]])
  }

  # Validate explicit options
  invalid <- setdiff(retain, RETAIN_OPTIONS)
  if (length(invalid) > 0) {
    cli::cli_warn("Unknown retention options ignored: {invalid}")
  }

  intersect(retain, RETAIN_OPTIONS)
}

## Fit Result Retention --------------------------------------------------------

#' Apply retention policy to fit result
#'
#' Removes fields from fit result based on retention policy to reduce
#' memory footprint. This function modifies the fit_result in place
#' by setting unwanted fields to NULL.
#'
#' @param fit_result A bayesim_fit_result object
#' @param retain Character vector of retention options specifying what to keep
#'
#' @return Modified bayesim_fit_result object with non-retained fields removed
#'
#' @export
apply_fit_retention <- function(fit_result, retain) {
  if (!"fit" %in% retain) {
    fit_result$fit <- NULL
  }
  if (!"draws" %in% retain) {
    fit_result$draws <- NULL
  }
  if (!"diagnostics" %in% retain) {
    fit_result$diagnostics <- NULL
  }
  if (!"warnings" %in% retain) {
    fit_result$warnings <- character()
  }

  fit_result
}

## Task Result Retention -------------------------------------------------------

#' Apply retention policy to task result
#'
#' Adds optional retained fields from the fit result and data bundle
#' to the task result based on the retention policy.
#'
#' @param task_result A bayesim_task_result object to modify
#' @param fit_result A bayesim_fit_result object (before retention applied)
#' @param data_bundle Data bundle from generator containing train/test data
#' @param retain Character vector of retention options specifying what to keep
#'
#' @return Modified bayesim_task_result object with retained fields added
#'
#' @export
apply_task_retention <- function(task_result, fit_result, data_bundle, retain) {
  # Add optional retained fields
  if ("draws" %in% retain && !is.null(fit_result$draws)) {
    task_result$draws <- fit_result$draws
  }
  if ("fit" %in% retain && !is.null(fit_result$fit)) {
    task_result$fit <- fit_result$fit
  }
  if ("data" %in% retain) {
    task_result$data <- list(
      train = data_bundle$train,
      test = data_bundle$test
    )
  }

  task_result
}

## Size Estimation -------------------------------------------------------------

#' Estimate object size
#'
#' Estimates the memory size of an R object in bytes.
#'
#' @param x An R object to estimate
#'
#' @return Size in bytes as numeric
#'
#' @export
#'
#' @examples
#' estimate_size(1:1000)
#' estimate_size(iris)
estimate_size <- function(x) {
  as.numeric(utils::object.size(x))
}

#' Check if result row exceeds size threshold
#'
#' Determines whether a task result exceeds a specified memory threshold,
#' useful for deciding whether to externalize large artifacts.
#'
#' @param task_result A bayesim_task_result object
#' @param threshold_bytes Maximum allowed size in bytes. Default is 5 MB.
#'
#' @return TRUE if the task_result exceeds the threshold, FALSE otherwise
#'
#' @export
#'
#' @examples
#' result <- list(metrics = list(rmse = 0.1), diagnostics = list(rhat = 1.01))
#' exceeds_size_threshold(result)
#' exceeds_size_threshold(result, threshold_bytes = 10)
exceeds_size_threshold <- function(
  task_result,
  threshold_bytes = 5 * 1024 * 1024
) {
  estimate_size(task_result) > threshold_bytes
}

## Externalization -------------------------------------------------------------

#' Externalize large artifact
#'
#' Moves a large artifact to an external RDS file and replaces the in-memory
#' object with a pointer containing metadata about the externalized file.
#'
#' @param artifact An R object to externalize
#' @param artifacts_dir Directory path where artifacts should be stored
#' @param task_id Task identifier used to construct the filename
#' @param field_name Name of the field being externalized (e.g., "draws", "fit")
#'
#' @return A list with pointer information containing:
#'   \itemize{
#'     \item external: TRUE (indicating this is an external pointer)
#'     \item path: Full path to the external file
#'     \item hash: Hash of the artifact for integrity checking
#'     \item size: Size of the artifact in bytes
#'   }
#'
#' @export
#'
#' @seealso [write_rds_atomic()], [compute_hash()]
externalize_artifact <- function(artifact, artifacts_dir, task_id, field_name) {
  dir.create(artifacts_dir, recursive = TRUE, showWarnings = FALSE)

  filename <- paste0(task_id, "__", field_name, ".rds")
  filepath <- file.path(artifacts_dir, filename)

  write_rds_atomic(artifact, filepath)

  list(
    external = TRUE,
    path = filepath,
    hash = compute_hash(artifact),
    size = estimate_size(artifact)
  )
}
