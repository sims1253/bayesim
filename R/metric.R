#' @title Metric Abstract Class
#' @description Abstract base class for defining simulation metrics in bayesim.
#'   Metrics compute summary statistics or diagnostic values from fitted model
#'   results and data bundles. `Metric` is abstract: it has no direct
#'   instances, so calling `Metric()` errors. Subclasses are created with
#'   `S7::new_class(..., parent = Metric)` and instantiated directly. Because
#'   the parent is abstract, S7 honors the defaults a subclass declares for
#'   the inherited `name`/`needs`/`required`/`summary_type`/`schema`
#'   properties when the subclass constructor is called.
#'
#' @param name Character identifier for the metric. Used as a prefix when
#'   flattening metric output to column names.
#' @param needs Character vector of required capabilities from the fitter.
#'   Common values include "predictions", "log_lik", "loo", "epred". The metric
#'   will only receive these values in the context if the fitter provides
#'   them. `epred` is delivered as `context$loo_epred` inside the LOO context,
#'   so declare it alongside `"loo"` (as `rmse_loo_metric()` does); a metric
#'   declaring `needs = "epred"` alone never receives the matrix.
#' @param required Logical indicating whether metric failure causes task
#'   failure. If TRUE, an error in computing this metric will propagate and
#'   fail the entire task. If FALSE (default), metric failure results in
#'   NA values being recorded.
#' @param summary_type Character; how [summarize_simulation()] aggregates this
#'   metric's flattened columns: `"mean"` (default, sd/sqrt(n) MCSE),
#'   `"proportion"` (coverage-style sqrt(p(1-p)/n) MCSE), or `"none"`.
#' @param schema Named list of field-level metadata. Each emitted field can
#'   declare a `role` (`estimate`, `binary`, `count`, `diagnostic`, `rank`, or
#'   `artifact`), an `aggregation` (`mean`, `proportion`, or `none`), an MCSE
#'   method (`sd`, `binomial`, or `none`), and optional `nominal`, `units`, or
#'   `dimension` metadata. `summary_type` remains supported as a compatibility
#'   default for metrics that do not declare a schema.
#'
#' @return An S7 class object representing the abstract Metric base class.
#'   Construct subclasses directly (e.g. `MyMetric()`); do not call `Metric()`.
#'
#' @section Methods:
#' The `compute_metric()` S7 generic must be implemented by subclasses.
#' \describe{
#'   \item{compute_metric(metric, fit_result, data_bundle, context, task_ctx)}{
#'     Compute metric values from a fitted model result. This method must be
#'     implemented by subclasses.
#'
#'     \itemize{
#'       \item metric: The Metric S7 object
#'       \item fit_result: A bayesim_fit_result object containing the fitted
#'         model output (draws, diagnostics, etc.)
#'       \item data_bundle: A list containing data-related objects including
#'         train (training data), test (test data if applicable), response
#'         (response variable), true_params (true parameter values if known)
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
#' S7::method(compute_metric, RMSEMetric) <- function(
#'   metric, fit_result, data_bundle, context, task_ctx
#' ) {
#'   preds <- context$predictions$predicted_mean
#'   actual <- data_bundle$test[[data_bundle$response]]
#'   list(
#'     value = sqrt(mean((preds - actual)^2)),
#'     n_obs = length(actual)
#'   )
#' }
#'
#' # Construct the subclass directly; the declared defaults are honored.
#' RMSEMetric()@name     # "rmse"
#' RMSEMetric()@needs    # "predictions"
#' RMSEMetric()@required # FALSE
#'
#' @export
Metric <- S7::new_class(
  "Metric",
  abstract = TRUE,
  properties = list(
    # Subclasses provide the concrete default (Metric is abstract, so S7
    # honors subclass-declared defaults at construction). The public
    # validator checks the value after subclass construction, avoiding
    # premature parent validation.
    name = S7::new_property(S7::class_character),
    needs = S7::new_property(S7::class_character, default = character()),
    required = S7::new_property(S7::class_logical, default = FALSE),
    # E4: how summarize_simulation aggregates this metric's flattened columns.
    # "mean" (default) — mean/sd/sqrt(n) MCSE; "proportion" — coverage-style
    # sqrt(p(1-p)/n) MCSE; "none" — do not aggregate (e.g. per-task ranks).
    summary_type = S7::new_property(
      S7::class_character,
      default = "mean",
      validator = function(value) validate_metric_summary_type(value)
    ),
    schema = S7::new_property(
      S7::class_list,
      default = list(),
      validator = function(value) validate_metric_schema(value)
    )
  )
)

# =============================================================================
# Metric metadata and shared validation
# =============================================================================

METRIC_SCHEMA_ROLES <- c(
  "estimate",
  "binary",
  "count",
  "diagnostic",
  "rank",
  "artifact"
)
METRIC_SCHEMA_AGGREGATIONS <- c("mean", "proportion", "none")
METRIC_SCHEMA_MCSE <- c("sd", "binomial", "none")

# S7 validators return NULL for valid values and a short message otherwise.
# Keeping these validators here makes built-in metric constructors reject bad
# configuration before a run starts, while validate_metric() provides the same
# diagnostics for package-external S7 subclasses.
validate_metric_name <- function(value) {
  if (
    !is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !nzchar(value)
  ) {
    return("name must be a non-empty character scalar")
  }
  if (grepl("__", value, fixed = TRUE)) {
    return("name must not contain the '__' flattening separator")
  }
  NULL
}

validate_metric_summary_type <- function(value) {
  if (
    !is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !value %in% METRIC_SCHEMA_AGGREGATIONS
  ) {
    paste0(
      "summary_type must be one of: ",
      paste(METRIC_SCHEMA_AGGREGATIONS, collapse = ", ")
    )
  }
}

validate_interval_probability <- function(value, label = "probability") {
  if (
    !is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value) ||
      value <= 0 ||
      value >= 1
  ) {
    paste0(label, " must be a single finite numeric value in (0, 1)")
  }
}

validate_metric_schema <- function(schema) {
  if (is.null(schema)) {
    return(NULL)
  }
  if (!is.list(schema)) {
    return("schema must be a named list of field metadata")
  }
  if (length(schema) == 0L) {
    return(NULL)
  }
  nms <- names(schema)
  if (
    is.null(nms) ||
      anyNA(nms) ||
      !all(nzchar(nms)) ||
      anyDuplicated(nms) > 0L
  ) {
    return("schema must have unique, non-empty field names")
  }
  for (field in nms) {
    meta <- schema[[field]]
    if (is.character(meta) && length(meta) == 1L) {
      if (!meta %in% METRIC_SCHEMA_AGGREGATIONS) {
        return(paste0("schema field '", field, "' has an invalid aggregation"))
      }
      next
    }
    if (!is.list(meta)) {
      return(paste0("schema field '", field, "' must be a metadata list"))
    }
    unknown <- setdiff(
      names(meta),
      c(
        "role",
        "aggregation",
        "mcse",
        "nominal",
        "dimension",
        "units",
        "externalize"
      )
    )
    if (length(unknown) > 0L) {
      return(paste0(
        "schema field '",
        field,
        "' has unknown metadata: ",
        paste(unknown, collapse = ", ")
      ))
    }
    if (
      !is.null(meta$role) &&
        (!is.character(meta$role) ||
          length(meta$role) != 1L ||
          is.na(meta$role) ||
          !meta$role %in% METRIC_SCHEMA_ROLES)
    ) {
      return(paste0("schema field '", field, "' has an invalid role"))
    }
    if (
      !is.null(meta$aggregation) &&
        (!is.character(meta$aggregation) ||
          length(meta$aggregation) != 1L ||
          is.na(meta$aggregation) ||
          !meta$aggregation %in% METRIC_SCHEMA_AGGREGATIONS)
    ) {
      return(paste0("schema field '", field, "' has an invalid aggregation"))
    }
    if (
      !is.null(meta$mcse) &&
        (!is.character(meta$mcse) ||
          length(meta$mcse) != 1L ||
          is.na(meta$mcse) ||
          !meta$mcse %in% METRIC_SCHEMA_MCSE)
    ) {
      return(paste0("schema field '", field, "' has an invalid mcse method"))
    }
    if (!is.null(meta$nominal)) {
      probability_error <- validate_interval_probability(
        meta$nominal,
        "nominal"
      )
      if (!is.null(probability_error)) {
        return(paste0("schema field '", field, "': ", probability_error))
      }
    }
    if (
      !is.null(meta$dimension) &&
        (!is.character(meta$dimension) ||
          length(meta$dimension) != 1L ||
          is.na(meta$dimension) ||
          !nzchar(meta$dimension))
    ) {
      return(paste0("schema field '", field, "' has an invalid dimension"))
    }
    if (
      !is.null(meta$units) &&
        (!is.character(meta$units) ||
          length(meta$units) != 1L ||
          is.na(meta$units))
    ) {
      return(paste0("schema field '", field, "' has invalid units"))
    }
    if (
      !is.null(meta$externalize) &&
        (!is.logical(meta$externalize) ||
          length(meta$externalize) != 1L ||
          is.na(meta$externalize))
    ) {
      return(paste0("schema field '", field, "' has invalid externalize flag"))
    }
  }
  NULL
}

# Return normalized metadata for a field. Character shorthand means an
# aggregation rule. A caller can pass a Metric object from any subclass; old
# subclasses that only declare summary_type continue to work.
metric_field_metadata <- function(metric, field = NULL) {
  schema <- tryCatch(metric@schema, error = function(e) list())
  summary_type <- tryCatch(metric@summary_type, error = function(e) "mean")
  if (!is.null(field) && length(schema) > 0L && field %in% names(schema)) {
    meta <- schema[[field]]
    if (is.character(meta)) meta <- list(aggregation = meta)
  } else {
    meta <- list()
  }
  if (is.null(meta$aggregation)) {
    meta$aggregation <- summary_type
  }
  if (is.null(meta$mcse)) {
    meta$mcse <- if (identical(meta$aggregation, "proportion")) {
      "binomial"
    } else if (identical(meta$aggregation, "none")) {
      "none"
    } else {
      "sd"
    }
  }
  meta
}

# Prediction metrics must never rely on R's vector recycling. A malformed
# fitter result is a contract violation, not a shorter prediction problem;
# fail at the metric seam so the worker can record a clear metric failure.
validate_prediction_vectors <- function(actual, predicted, metric_name) {
  is_vector <- function(x) is.atomic(x) && is.null(dim(x))
  if (!is_vector(actual) || !is.numeric(actual)) {
    stop(bayesim_metric_error(
      sprintf("Metric '%s' requires a numeric response vector", metric_name)
    ))
  }
  if (!is_vector(predicted) || !is.numeric(predicted)) {
    stop(bayesim_metric_error(
      sprintf("Metric '%s' requires a numeric prediction vector", metric_name)
    ))
  }
  if (length(actual) == 0L || length(predicted) == 0L) {
    stop(bayesim_metric_error(
      sprintf(
        "Metric '%s' requires non-empty response and prediction vectors",
        metric_name
      )
    ))
  }
  if (length(actual) != length(predicted)) {
    stop(bayesim_metric_error(
      sprintf(
        "Metric '%s' received %d responses but %d predictions; lengths must match",
        metric_name,
        length(actual),
        length(predicted)
      )
    ))
  }
  invisible(NULL)
}

# =============================================================================
# S7 Generic for Metric compute method
# =============================================================================

#' @title Compute Metric Values
#' @description
#' Compute metric values from a fitted model result. This generic must be
#' implemented by Metric subclasses. Named `compute_metric` (rather than
#' `compute`) so that bayesim does not mask [dplyr::compute].
#'
#' @param metric An S7 Metric object
#' @param fit_result A bayesim_fit_result object containing the fitted model output
#' @param data_bundle A list containing data-related objects
#' @param context A list with precomputed values based on the metric's `needs`
#' @param task_ctx A list with task identification information
#'
#' @return A named list with metric values conforming to the metric output schema
#' @export
compute_metric <- S7::new_generic(
  "compute_metric",
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
#' @return Invisible `output` if validation passes. Otherwise, an error is raised.
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
#' \dontrun{
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
#' }
#' @keywords internal
validate_metric_output <- function(output, metric_name) {
  if (!is.list(output)) {
    stop(
      bayesim_metric_error(
        sprintf(
          "Metric '%s' output must be a list, got %s",
          metric_name,
          class(output)[1]
        )
      )
    )
  }

  if (length(output) == 0) {
    stop(
      bayesim_metric_error(
        sprintf("Metric '%s' output cannot be empty", metric_name)
      )
    )
  }

  nms <- names(output)
  if (is.null(nms) || any(nms == "" | is.na(nms))) {
    unnamed_idx <- which(is.null(nms) | nms == "" | is.na(nms))
    stop(
      bayesim_metric_error(
        sprintf(
          "Metric '%s' output has unnamed or empty-named elements at positions: %s",
          metric_name,
          paste(unnamed_idx, collapse = ", ")
        )
      )
    )
  }
  if (anyDuplicated(nms) > 0L) {
    stop(
      bayesim_metric_error(
        sprintf(
          "Metric '%s' output has duplicate field names: %s",
          metric_name,
          paste(unique(nms[duplicated(nms)]), collapse = ", ")
        )
      )
    )
  }

  for (nm in nms) {
    val <- output[[nm]]

    if (is.null(val)) {
      stop(
        bayesim_metric_error(
          sprintf(
            "Metric '%s' output element '%s' is NULL (not allowed)",
            metric_name,
            nm
          )
        )
      )
    }

    is_scalar_atomic <- (is.logical(val) ||
      is.integer(val) ||
      is.double(val) ||
      is.character(val)) &&
      length(val) == 1

    is_named_numeric_vector <- is.double(val) &&
      length(val) > 1 &&
      !is.null(names(val)) &&
      all(names(val) != "" & !is.na(names(val))) &&
      anyDuplicated(names(val)) == 0L

    if (!is_scalar_atomic && !is_named_numeric_vector) {
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
        bayesim_metric_error(
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
          )
        )
      )
    }
  }

  invisible(output)
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
#' \dontrun{
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
#' }
#' @keywords internal
flatten_metric_output <- function(output, metric_name) {
  validate_metric_output(output, metric_name)

  result <- list()

  for (nm in names(output)) {
    val <- output[[nm]]

    # A named numeric vector is expanded to <metric>__<field>__<subname> per
    # element, EVEN when length 1 — so single-parameter results (e.g. a
    # one-variable rank/coverage by_param) still carry the <param> suffix and
    # downstream consumers grepping for <metric>__<field>__<param> find them.
    # Bare scalars (no names) collapse to <metric>__<field>. validate_metric_
    # output() only admits these two shapes: scalar atomics (length 1) and
    # named numeric vectors (length >= 1), so no other branch is needed.
    if (!is.null(names(val)) && is.numeric(val) && length(val) >= 1L) {
      for (sub_nm in names(val)) {
        result[[paste0(metric_name, "__", nm, "__", sub_nm)]] <- val[[sub_nm]]
      }
    } else {
      result[[paste0(metric_name, "__", nm)]] <- val
    }
  }

  result
}
