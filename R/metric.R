#' @title Metric Abstract Class
#' @description Abstract base class for defining simulation metrics in bayesim.
#'   Metrics compute summary statistics or diagnostic values from fitted model
#'   results and data bundles.
#'
#' @name Metric-class
#' @return An S7 class object representing a Metric.
#'
#' @section Properties:
#' \describe{
#'   \item{name}{Character identifier for the metric. Used as a prefix when
#'     flattening metric output to column names.}
#'   \item{needs}{Character vector of required capabilities from the fitter.
#'     Common values include "predictions", "log_lik", "loo". The metric
#'     will only receive these values in the context if the fitter provides
#'     them.}
#'   \item{required}{Logical indicating whether metric failure causes task
#'     failure. If TRUE, an error in computing this metric will propagate and
#'     fail the entire task. If FALSE (default), metric failure results in
#'     NA values being recorded.}
#' }
#'
#' @section Methods:
#' The `compute()` S7 generic must be implemented by subclasses.
#' \describe{
#'   \item{compute(metric, fit_result, data_bundle, context, task_ctx)}{
#'     Compute metric values from a fitted model result. This method must be
#'     implemented by subclasses.
#'
#'     \itemize{
#'       \item metric: The Metric S7 object
#'       \item fit_result: A bayesim_fit_result object containing the fitted
#'         model output (draws, diagnostics, etc.)
#'       \item data_bundle: A list containing data-related objects including
#'         train (training data), test (test data if applicable), response
#'         (response variable), true_params (true parameter values if known),
#'         references (reference values), etc.
#'       \item context: A list with precomputed values based on the metric's
#'         `needs` property. May include predictions, log_lik (log-likelihood
#'         values), loo (leave-one-out cross-validation results), etc.
#'       \item task_ctx: A list with task identification information including
#'         task_id, data_idx, fit_idx, rep_idx for tracking and debugging.
#'     }
#'
#'     Returns: A named list with metric values. Names must be non-empty strings.
#'     Values must be one of:
#'     \itemize{
#'       \item scalar atomic (logical, integer, double, character)
#'       \item named numeric vector
#'     }
#'
#'     No nested data frames or matrices are allowed in the output.
#'   }
#' }
#'
#' @section Metric Output Schema:
#' The compute method must return a named list conforming to the following
#' schema:
#' \itemize{
#'   \item All elements must have non-empty names
#'   \item Values must be scalar atomic types or named numeric vectors
#'   \item No nested data frames or matrices allowed
#'   \item The engine flattens output with prefix `<metric_name>__<field>`
#' }
#'
#' @examples
#' # Define a custom RMSE metric
#' RMSEMetric <- S7::new_class(
#'   "RMSEMetric",
#'   parent = Metric,
#'   properties = list(
#'     name = S7::new_property(S7::class_character, default = "rmse"),
#'     needs = S7::new_property(S7::class_character, default = "predictions"),
#'     required = S7::new_property(S7::class_logical, default = FALSE)
#'   )
#' )
#'
#' S7::method(compute, RMSEMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
#'   preds <- context$predictions
#'   actual <- data_bundle$test[[data_bundle$response]]
#'   list(
#'     value = sqrt(mean((preds - actual)^2)),
#'     n_obs = length(actual)
#'   )
#' }
#'
#' @export
Metric <- S7::new_class(
  "Metric",
  properties = list(
    name = S7::new_property(S7::class_character),
    needs = S7::new_property(S7::class_character, default = character()),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

# =============================================================================
# S7 Generic for Metric compute method
# =============================================================================

#' @title Compute Metric Values
#' @description
#' Compute metric values from a fitted model result. This generic must be
#' implemented by Metric subclasses.
#'
#' @param metric An S7 Metric object
#' @param fit_result A bayesim_fit_result object containing the fitted model output
#' @param data_bundle A list containing data-related objects
#' @param context A list with precomputed values based on the metric's `needs`
#' @param task_ctx A list with task identification information
#'
#' @return A named list with metric values conforming to the metric output schema
#' @export
compute <- S7::new_generic(
  "compute",
  "metric",
  function(metric, fit_result, data_bundle, context, task_ctx) {
    S7::S7_dispatch()
  }
)


#' @title Validate Metric Output
#' @description Validates that metric output conforms to the expected schema.
#'   Metric output must be a named list where all names are non-empty strings
#'   and values are either scalar atomic types or named numeric vectors.
#'
#' @param output The output list to validate.
#' @param metric_name Character string identifying the metric (for error messages).
#'
#' @return Invisible NULL if validation passes. Otherwise, an error is raised.
#'
#' @section Validation Rules:
#' \itemize{
#'   \item output must be a list
#'   \item all elements must have non-empty names
#'   \item values must be one of:
#'     \itemize{
#'       \item scalar logical, integer, double, or character
#'       \item named numeric vector
#'     }
#'   \item no NULL values allowed
#'   \item no nested lists, data frames, or matrices allowed
#' }
#'
#' @examples
#' # Valid output
#' validate_metric_output(list(rmse = 0.5, n_obs = 100L), "rmse")
#'
#' # Valid output with named vector
#' validate_metric_output(
#'   list(params = c(alpha = 0.1, beta = 0.2)),
#'   "param_estimates"
#' )
#'
#' \dontrun{
#' # Invalid: unnamed element
#' validate_metric_output(list(0.5), "rmse")  # Error
#'
#' # Invalid: nested list
#' validate_metric_output(list(nested = list(a = 1)), "rmse")  # Error
#' }
#'
#' @export
validate_metric_output <- function(output, metric_name) {
  # Check that output is a list
  if (!is.list(output)) {
    stop(
      sprintf(
        "Metric '%s' output must be a list, got %s",
        metric_name,
        class(output)[1]
      ),
      call. = FALSE
    )
  }

  # Check for empty output
  if (length(output) == 0) {
    stop(
      sprintf("Metric '%s' output cannot be empty", metric_name),
      call. = FALSE
    )
  }

  # Get names and check they are non-empty
  nms <- names(output)
  if (is.null(nms) || any(nms == "" | is.na(nms))) {
    unnamed_idx <- which(is.null(nms) | nms == "" | is.na(nms))
    stop(
      sprintf(
        "Metric '%s' output has unnamed or empty-named elements at positions: %s",
        metric_name,
        paste(unnamed_idx, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Check each value
  for (nm in nms) {
    val <- output[[nm]]

    # Check for NULL
    if (is.null(val)) {
      stop(
        sprintf(
          "Metric '%s' output element '%s' is NULL (not allowed)",
          metric_name,
          nm
        ),
        call. = FALSE
      )
    }

    # Check for allowed types: scalar atomic or named numeric vector
    is_scalar_atomic <- (is.logical(val) ||
      is.integer(val) ||
      is.double(val) ||
      is.character(val)) &&
      length(val) == 1

    is_named_numeric_vector <- is.double(val) &&
      length(val) > 1 &&
      !is.null(names(val)) &&
      all(names(val) != "" & !is.na(names(val)))

    if (!is_scalar_atomic && !is_named_numeric_vector) {
      # Determine what we got for error message
      if (is.list(val)) {
        type_desc <- "list"
      } else if (is.data.frame(val)) {
        type_desc <- "data.frame"
      } else if (is.matrix(val)) {
        type_desc <- "matrix"
      } else if (is.factor(val)) {
        type_desc <- "factor"
      } else if (
        is.double(val) && (is.null(names(val)) || any(names(val) == ""))
      ) {
        type_desc <- "unnamed numeric vector"
      } else if (
        (is.logical(val) || is.integer(val) || is.character(val)) &&
          length(val) > 1
      ) {
        type_desc <- paste0("vector of type ", typeof(val))
      } else {
        type_desc <- paste(class(val), collapse = "/")
      }

      stop(
        sprintf(
          paste0(
            "Metric '%s' output element '%s' must be either:\n",
            "  - a scalar atomic (logical, integer, double, character), or\n",
            "  - a named numeric vector\n",
            "Got: %s with length %d"
          ),
          metric_name,
          nm,
          type_desc,
          length(val)
        ),
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}

#' @title Flatten Metric Output
#' @description Flattens a metric output list into a single-level named list
#'   with prefixed column names. This is used by the simulation engine to
#'   convert metric outputs into flat columns for results collection.
#'
#' @param output The metric output list to flatten.
#' @param metric_name Character string identifying the metric, used as prefix.
#'
#' @return A named list with flattened names. Scalar values get names like
#'   `<metric_name>__<field>`. Named numeric vectors get expanded to
#'   `<metric_name>__<field>__<subname>`.
#'
#' @section Flattening Rules:
#' \itemize{
#'   \item Scalar values become `<metric_name>__<field>`
#'   \item Named numeric vectors expand to `<metric_name>__<field>__<subname>`
#'     for each element
#'   \item The double underscore `__` is used as a separator consistently
#' }
#'
#' @examples
#' # Scalar values
#' flatten_metric_output(list(rmse = 0.5, n_obs = 100L), "my_metric")
#' # Returns: list(my_metric__rmse = 0.5, my_metric__n_obs = 100L)
#'
#' # Named numeric vectors are expanded
#' flatten_metric_output(
#'   list(params = c(alpha = 0.1, beta = 0.2, gamma = 0.3)),
#'   "estimates"
#' )
#' # Returns: list(
#' #   estimates__params__alpha = 0.1,
#' #   estimates__params__beta = 0.2,
#' #   estimates__params__gamma = 0.3
#' # )
#'
#' @export
flatten_metric_output <- function(output, metric_name) {
  # Validate first
  validate_metric_output(output, metric_name)

  # Initialize result
  result <- list()

  # Process each element
  for (nm in names(output)) {
    val <- output[[nm]]

    # Check if scalar or named vector
    if (length(val) == 1) {
      # Scalar: just prefix with metric name
      result[[paste0(metric_name, "__", nm)]] <- val
    } else {
      # Named numeric vector: expand with subnames
      for (sub_nm in names(val)) {
        result[[paste0(metric_name, "__", nm, "__", sub_nm)]] <- val[[sub_nm]]
      }
    }
  }

  result
}
