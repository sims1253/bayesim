# One-time warning helper: metric compute() runs once per task, so a repeated
# condition (missing test set, unsupported fitter) would otherwise warn
# thousands of times per run. Warn once per key per run; .reset_warn_once()
# clears the flags at the start of each run_simulation().
.warn_once_env <- new.env(parent = emptyenv())

# S6: clear all one-time-warning flags. Called at the start of run_simulation()
# so that each run warns on its own conditions even within one session.
.reset_warn_once <- function() {
  rm(list = ls(.warn_once_env, all.names = TRUE), envir = .warn_once_env)
  invisible(NULL)
}

.warn_once <- function(key, ..., .envir = parent.frame()) {
  if (is.null(.warn_once_env[[key]])) {
    .warn_once_env[[key]] <- TRUE
    # .envir: cli glue expressions in `...` must resolve in the caller.
    cli::cli_warn(c(...), .envir = .envir)
  }
  invisible(NULL)
}

# E2: prediction-error metrics refuse to silently fall back to the training
# set. Warn once per metric name per session, naming the fix (provide a test
# set). In-sample prediction error presented as "rmse" is a trap.
.warn_no_test <- function(metric_name) {
  .warn_once(
    metric_name,
    "{metric_name}: no test set in data_bundle; returning NA.",
    i = "Provide a test set via the data generator to measure predictive error."
  )
}

#' Resolve requested cleaned var names to actual draws-matrix columns
#'
#' Maps each name in `vars` to itself, then to `b_<name>`, in the set of
#' `draws_colnames`, mirroring the both-directions lookup `.extract_truth()`
#' already performs. This fixes the F2 mismatch where generators strip the
#' `b_` prefix (so `vars_of_interest` holds cleaned names like `c("x",
#' "Intercept")`) but brms draws matrices keep it (`c("b_x","b_Intercept",
#' "sigma")`).
#'
#' Errors with a `bayesim_config_error` (same condition class as
#' `.extract_truth()`) when a requested var is genuinely absent. A silent NA
#' would corrupt SBC ranks and credible intervals without diagnostics.
#'
#' @param vars Character vector of cleaned names (e.g. `c("x","Intercept",
#'   "sigma")`).
#' @param draws_colnames Character vector of available draws-matrix column
#'   names.
#'
#' @return A **named** character vector: names are the cleaned `vars` (use
#'   these for output field naming), values are the actual draws column to
#'   read. Empty input returns `character(0)` (names preserved).
#'
#' @keywords internal
resolve_draw_columns <- function(vars, draws_colnames) {
  if (length(vars) == 0L) {
    out <- character(0)
    names(out) <- character(0)
    return(out)
  }
  hits <- vars
  names(hits) <- vars
  missing <- character(0)
  for (i in seq_along(vars)) {
    v <- vars[[i]]
    if (v %in% draws_colnames) {
      next
    } else if (paste0("b_", v) %in% draws_colnames) {
      hits[[i]] <- paste0("b_", v)
    } else {
      missing <- c(missing, v)
    }
  }
  if (length(missing) > 0L) {
    stop(bayesim_config_error(
      "vars_of_interest not found in model draws: " %+%
        paste(missing, collapse = ", ") %+%
        ". Available: " %+%
        paste(utils::head(draws_colnames, 20), collapse = ", ")
    ))
  }
  hits
}

#' @title RMSE Metric
#' @description Root Mean Square Error between predictions and true values.
#' @param name Character string naming the metric. Defaults to "rmse".
#' @param needs Character vector of required capabilities. Defaults to "predictions".
#' @param required Logical indicating if metric failure causes task failure. Defaults to FALSE.
#' @keywords internal
RmseMetric <- S7::new_class(
  "RmseMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(
      S7::class_character,
      default = "rmse",
      validator = validate_metric_name
    ),
    needs = S7::new_property(S7::class_character, default = "predictions"),
    required = S7::new_property(S7::class_logical, default = FALSE),
    schema = S7::new_property(
      S7::class_list,
      default = list(
        value = list(role = "estimate", aggregation = "mean", mcse = "sd"),
        n_obs = list(role = "count", aggregation = "none", mcse = "none")
      ),
      validator = function(value) validate_metric_schema(value)
    )
  )
)

#' @rdname RmseMetric
#' @description Constructor for RmseMetric.
#' @param name Character string naming the metric. Defaults to "rmse".
#' @return An RmseMetric object.
#' @export
#' @examples
#' pred_rmse_metric()
pred_rmse_metric <- function(name = "rmse") {
  RmseMetric(
    name = name,
    needs = "predictions",
    required = FALSE,
    schema = list(
      value = list(role = "estimate", aggregation = "mean", mcse = "sd"),
      n_obs = list(role = "count", aggregation = "none", mcse = "none")
    )
  )
}

S7::method(compute_metric, RmseMetric) <- function(
  metric,
  fit_result,
  data_bundle,
  context,
  task_ctx
) {
  if (is.null(context$predictions)) {
    return(list(value = NA_real_, n_obs = NA_integer_))
  }
  # E2: prediction-error metrics must NOT silently fall back to the training
  # set — in-sample error presented as "rmse" is a trap for this audience.
  if (is.null(data_bundle$test)) {
    .warn_no_test("pred_rmse_metric")
    return(list(value = NA_real_, n_obs = NA_integer_))
  }

  test_data <- data_bundle$test
  actual <- test_data[[data_bundle$response]]
  predicted <- context$predictions$predicted_mean
  validate_prediction_vectors(actual, predicted, metric@name)

  list(
    value = sqrt(mean((predicted - actual)^2)),
    n_obs = length(actual)
  )
}

#' @title Bias Metric
#' @description Mean bias of predictions.
#'
#' @param name Character string naming the metric. Defaults to "bias".
#' @param needs Character vector of required capabilities. Defaults to "predictions".
#' @param required Logical indicating if metric failure causes task failure. Defaults to FALSE.
#'
#' @return A BiasMetric object.
#' @keywords internal
BiasMetric <- S7::new_class(
  "BiasMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(
      S7::class_character,
      default = "bias",
      validator = validate_metric_name
    ),
    needs = S7::new_property(S7::class_character, default = "predictions"),
    required = S7::new_property(S7::class_logical, default = FALSE),
    schema = S7::new_property(
      S7::class_list,
      default = list(
        value = list(role = "estimate", aggregation = "mean", mcse = "sd")
      ),
      validator = function(value) validate_metric_schema(value)
    )
  )
)

#' @rdname BiasMetric
#' @description Constructor for BiasMetric.
#' @param name Character string naming the metric. Defaults to "bias".
#' @return A BiasMetric object.
#' @export
pred_bias_metric <- function(name = "bias") {
  BiasMetric(
    name = name,
    needs = "predictions",
    required = FALSE,
    schema = list(
      value = list(role = "estimate", aggregation = "mean", mcse = "sd")
    )
  )
}

S7::method(compute_metric, BiasMetric) <- function(
  metric,
  fit_result,
  data_bundle,
  context,
  task_ctx
) {
  if (is.null(context$predictions)) {
    return(list(value = NA_real_))
  }
  # E2: no silent training-set fallback.
  if (is.null(data_bundle$test)) {
    .warn_no_test("pred_bias_metric")
    return(list(value = NA_real_))
  }

  test_data <- data_bundle$test
  actual <- test_data[[data_bundle$response]]
  predicted <- context$predictions$predicted_mean
  validate_prediction_vectors(actual, predicted, metric@name)

  list(
    value = mean(predicted - actual)
  )
}

#' @title Coverage Metric
#' @description Coverage of true parameter values within credible intervals.
#' @param name Character string naming the metric. Defaults to "coverage".
#' @param needs Character vector of required capabilities. Defaults to empty character vector.
#' @param required Logical indicating if metric failure causes task failure. Defaults to FALSE.
#' @param prob Numeric probability for the credible interval width. Defaults to 0.95.
#' @keywords internal
CoverageMetric <- S7::new_class(
  "CoverageMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(
      S7::class_character,
      default = "coverage",
      validator = validate_metric_name
    ),
    needs = S7::new_property(S7::class_character, default = character()),
    required = S7::new_property(S7::class_logical, default = FALSE),
    # E4: coverage columns are proportions -> sqrt(p(1-p)/n) MCSE.
    summary_type = S7::new_property(
      S7::class_character,
      default = "proportion",
      validator = validate_metric_summary_type
    ),
    prob = S7::new_property(
      S7::class_numeric,
      default = 0.95,
      validator = function(value) validate_interval_probability(value, "prob")
    ),
    schema = S7::new_property(
      S7::class_list,
      default = list(
        mean = list(role = "estimate", aggregation = "mean", mcse = "sd"),
        by_param = list(
          role = "binary",
          aggregation = "proportion",
          mcse = "binomial"
        )
      ),
      validator = function(value) validate_metric_schema(value)
    )
  )
)

#' @rdname CoverageMetric
#' @description Constructor for CoverageMetric.
#' @param name Character string naming the metric. Defaults to "coverage".
#' @param prob Numeric probability for the credible interval width. Defaults to 0.95.
#' @return A CoverageMetric object.
#' @export
coverage_metric <- function(name = "coverage", prob = 0.95) {
  tryCatch(
    CoverageMetric(
      name = name,
      needs = character(),
      required = FALSE,
      summary_type = "proportion",
      prob = prob,
      schema = list(
        mean = list(role = "estimate", aggregation = "mean", mcse = "sd"),
        by_param = list(
          role = "binary",
          aggregation = "proportion",
          mcse = "binomial",
          nominal = prob
        )
      )
    ),
    error = function(e) {
      stop(bayesim_config_error(conditionMessage(e)))
    }
  )
}

S7::method(compute_metric, CoverageMetric) <- function(
  metric,
  fit_result,
  data_bundle,
  context,
  task_ctx
) {
  if (is.null(fit_result$draws) || is.null(data_bundle$true_params)) {
    return(list(value = NA_real_))
  }

  draws <- fit_result$draws
  true_params <- data_bundle$true_params
  vars_of_interest <- data_bundle$vars_of_interest

  lower_q <- (1 - metric@prob) / 2
  upper_q <- 1 - lower_q

  mapped <- resolve_draw_columns(vars_of_interest, colnames(draws))
  coverage <- vapply(
    names(mapped),
    function(vname) {
      if (!(vname %in% names(true_params))) {
        return(NA_real_)
      }
      col <- mapped[[vname]]
      ci <- quantile(draws[, col], c(lower_q, upper_q))
      as.numeric(true_params[[vname]] >= ci[1] && true_params[[vname]] <= ci[2])
    },
    FUN.VALUE = numeric(1),
    USE.NAMES = TRUE
  )

  list(
    mean = mean(coverage, na.rm = TRUE),
    by_param = coverage
  )
}
