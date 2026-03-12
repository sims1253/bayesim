# Global metric registry environment
.metric_registry <- new.env(parent = emptyenv())

#' Register a metric
#'
#' Adds a metric to the global registry for lookup by name.
#'
#' @param metric An S7 Metric object
#' @param overwrite Logical; if TRUE, overwrites existing metric with same name
#'
#' @return Invisible NULL
#'
#' @noRd
register_metric <- function(metric, overwrite = FALSE) {
  if (!S7::S7_inherits(metric)) {
    cli::cli_abort("metric must be an S7 object")
  }

  # Check it's actually a Metric subclass
  if (!S7::S7_inherits(metric, Metric)) {
    cli::cli_abort("metric must inherit from the Metric S7 class")
  }

  name <- metric@name
  if (is.null(name) || name == "") {
    cli::cli_abort("metric must have a non-empty name")
  }

  if (!overwrite && !is.null(.metric_registry[[name]])) {
    cli::cli_abort(
      "Metric '{name}' is already registered. Use overwrite = TRUE to replace."
    )
  }

  .metric_registry[[name]] <- metric
  invisible(NULL)
}

#' Get a metric by name
#'
#' Retrieves a metric from the registry.
#'
#' @param name Character; the metric name
#'
#' @return The S7 Metric object, or NULL if not found
#'
#' @noRd
get_metric <- function(name) {
  if (!is.character(name) || length(name) != 1) {
    cli::cli_abort("name must be a single character string")
  }

  .metric_registry[[name]]
}

#' List all registered metrics
#'
#' @return Character vector of registered metric names
#'
#' @noRd
list_metrics <- function() {
  ls(.metric_registry)
}

#' Unregister a metric
#'
#' Removes a metric from the registry.
#'
#' @param name Character; the metric name to remove
#'
#' @return Invisible NULL
#'
#' @noRd
unregister_metric <- function(name) {
  if (!is.character(name) || length(name) != 1) {
    cli::cli_abort("name must be a single character string")
  }

  if (is.null(.metric_registry[[name]])) {
    cli::cli_warn("Metric '{name}' is not registered")
    return(invisible(NULL))
  }

  rm(list = name, envir = .metric_registry)
  invisible(NULL)
}

#' Resolve metrics to Metric objects
#'
#' Converts character metric names to Metric objects via the registry.
#'
#' @param metrics Character vector of names, or list of Metric objects
#'
#' @return List of Metric objects
#'
#' @keywords internal
resolve_metrics_from_registry <- function(metrics) {
  if (is.null(metrics)) {
    return(list())
  }

  if (is.character(metrics)) {
    resolved <- lapply(metrics, function(name) {
      m <- get_metric(name)
      if (is.null(m)) {
        cli::cli_abort("Metric '{name}' not found in registry")
      }
      m
    })
    return(resolved)
  }

  if (is.list(metrics)) {
    # Validate each is a Metric
    for (i in seq_along(metrics)) {
      m <- metrics[[i]]
      if (!S7::S7_inherits(m)) {
        cli::cli_abort("metrics[[{i}]] is not an S7 Metric object")
      }
    }
    return(metrics)
  }

  cli::cli_abort("metrics must be a character vector or list of Metric objects")
}

#' Clear all registered metrics
#'
#' Removes all metrics from the registry. Primarily for testing.
#'
#' @return Invisible NULL
#'
#' @keywords internal
clear_registry <- function() {
  rm(list = ls(.metric_registry), envir = .metric_registry)
  invisible(NULL)
}
