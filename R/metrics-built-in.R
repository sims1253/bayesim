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
      "vars_of_interest not found in model draws: "
        %+% paste(missing, collapse = ", ")
        %+% ". Available: "
        %+% paste(utils::head(draws_colnames, 20), collapse = ", ")
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
    name = S7::new_property(S7::class_character, default = "rmse"),
    needs = S7::new_property(S7::class_character, default = "predictions"),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname RmseMetric
#' @description Constructor for RmseMetric.
#' @param name Character string naming the metric. Defaults to "rmse".
#' @return An RmseMetric object.
#' @export
#' @examples
#' rmse_metric()
rmse_metric <- function(name = "rmse") {
  RmseMetric(
    name = name,
    needs = "predictions",
    required = FALSE
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
    return(list(value = NA_real_))
  }

  test_data <- data_bundle$test %||% data_bundle$train
  actual <- test_data[[data_bundle$response]]
  predicted <- context$predictions$predicted_mean

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
    name = S7::new_property(S7::class_character, default = "bias"),
    needs = S7::new_property(S7::class_character, default = "predictions"),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname BiasMetric
#' @description Constructor for BiasMetric.
#' @param name Character string naming the metric. Defaults to "bias".
#' @return A BiasMetric object.
#' @export
bias_metric <- function(name = "bias") {
  BiasMetric(
    name = name,
    needs = "predictions",
    required = FALSE
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

  test_data <- data_bundle$test %||% data_bundle$train
  actual <- test_data[[data_bundle$response]]
  predicted <- context$predictions$predicted_mean

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
    name = S7::new_property(S7::class_character, default = "coverage"),
    needs = S7::new_property(S7::class_character, default = character()),
    required = S7::new_property(S7::class_logical, default = FALSE),
    prob = S7::new_property(S7::class_numeric, default = 0.95)
  )
)

#' @rdname CoverageMetric
#' @description Constructor for CoverageMetric.
#' @param name Character string naming the metric. Defaults to "coverage".
#' @param prob Numeric probability for the credible interval width. Defaults to 0.95.
#' @return A CoverageMetric object.
#' @export
coverage_metric <- function(name = "coverage", prob = 0.95) {
  CoverageMetric(
    name = name,
    needs = character(),
    required = FALSE,
    prob = prob
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
  coverage <- vapply(names(mapped), function(vname) {
    if (!(vname %in% names(true_params))) {
      return(NA_real_)
    }
    col <- mapped[[vname]]
    ci <- quantile(draws[, col], c(lower_q, upper_q))
    as.numeric(true_params[[vname]] >= ci[1] && true_params[[vname]] <= ci[2])
  }, FUN.VALUE = numeric(1), USE.NAMES = TRUE)

  list(
    mean = mean(coverage, na.rm = TRUE),
    by_param = coverage
  )
}

#' @title Posterior Mean Metric
#' @description Posterior mean estimates for parameters.
#' @param name Character string naming the metric. Defaults to "posterior_mean".
#' @param needs Character vector of required capabilities. Defaults to empty character vector.
#' @param required Logical indicating if metric failure causes task failure. Defaults to FALSE.
#' @keywords internal
PosteriorMeanMetric <- S7::new_class(
  "PosteriorMeanMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "posterior_mean"),
    needs = S7::new_property(S7::class_character, default = character()),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname PosteriorMeanMetric
#' @description Constructor for PosteriorMeanMetric.
#' @param name Character string naming the metric. Defaults to "posterior_mean".
#' @return A PosteriorMeanMetric object.
#' @keywords internal
posterior_mean_metric <- function(name = "posterior_mean") {
  PosteriorMeanMetric(
    name = name,
    needs = character(),
    required = FALSE
  )
}

S7::method(compute_metric, PosteriorMeanMetric) <- function(
  metric,
  fit_result,
  data_bundle,
  context,
  task_ctx
) {
  if (is.null(fit_result$draws)) {
    return(list(value = NA_real_))
  }

  draws <- fit_result$draws
  voi <- data_bundle$vars_of_interest %||% colnames(draws)
  mapped <- resolve_draw_columns(voi, colnames(draws))

  means <- colMeans(draws[, mapped, drop = FALSE])
  names(means) <- names(mapped)

  as.list(means)
}
