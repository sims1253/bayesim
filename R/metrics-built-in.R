#' @title RMSE Metric
#' @description Root Mean Square Error between predictions and true values.
#' @param name Character string naming the metric. Defaults to "rmse".
#' @param needs Character vector of required capabilities. Defaults to "predictions".
#' @param required Logical indicating if metric failure causes task failure. Defaults to FALSE.
#' @export
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
#' @description `rmse_metric()` is a constructor function that creates an RmseMetric
#'   instance with appropriate defaults. Use this constructor to work around S7's
#'   property default inheritance issue.
#' @param name Character string naming the metric. Defaults to "rmse".
#' @return An RmseMetric object.
#' @export
#' @examples
#' rmse_metric()
#' rmse_metric(name = "my_rmse")
rmse_metric <- function(name = "rmse") {
  RmseMetric(
    name = name,
    needs = "predictions",
    required = FALSE
  )
}

S7::method(compute, RmseMetric) <- function(
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
#' @export
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
#' @description `bias_metric()` is a constructor function that creates a BiasMetric
#'   instance with appropriate defaults. Use this constructor to work around S7's
#'   property default inheritance issue.
#' @param name Character string naming the metric. Defaults to "bias".
#' @return A BiasMetric object.
#' @export
#' @examples
#' bias_metric()
#' bias_metric(name = "my_bias")
bias_metric <- function(name = "bias") {
  BiasMetric(
    name = name,
    needs = "predictions",
    required = FALSE
  )
}

S7::method(compute, BiasMetric) <- function(
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
#' @export
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
#' @description `coverage_metric()` is a constructor function that creates a CoverageMetric
#'   instance with appropriate defaults. Use this constructor to work around S7's
#'   property default inheritance issue.
#' @param name Character string naming the metric. Defaults to "coverage".
#' @param prob Numeric probability for the credible interval width. Defaults to 0.95.
#' @return A CoverageMetric object.
#' @export
#' @examples
#' coverage_metric()
#' coverage_metric(prob = 0.90)
coverage_metric <- function(name = "coverage", prob = 0.95) {
  CoverageMetric(
    name = name,
    needs = character(),
    required = FALSE,
    prob = prob
  )
}

S7::method(compute, CoverageMetric) <- function(
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

  coverage <- sapply(vars_of_interest, function(var) {
    if (!(var %in% colnames(draws)) || !(var %in% names(true_params))) {
      return(NA_real_)
    }

    ci <- quantile(draws[, var], c(lower_q, upper_q))
    as.numeric(true_params[var] >= ci[1] && true_params[var] <= ci[2])
  })

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
#' @export
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
#' @description `posterior_mean_metric()` is a constructor function that creates a
#'   PosteriorMeanMetric instance with appropriate defaults. Use this constructor to
#'   work around S7's property default inheritance issue.
#' @param name Character string naming the metric. Defaults to "posterior_mean".
#' @return A PosteriorMeanMetric object.
#' @export
#' @examples
#' posterior_mean_metric()
#' posterior_mean_metric(name = "my_posterior_mean")
posterior_mean_metric <- function(name = "posterior_mean") {
  PosteriorMeanMetric(
    name = name,
    needs = character(),
    required = FALSE
  )
}

S7::method(compute, PosteriorMeanMetric) <- function(
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
  vars_of_interest <- data_bundle$vars_of_interest %||% colnames(draws)

  means <- colMeans(draws[, vars_of_interest, drop = FALSE])

  as.list(means)
}

#' Register built-in metrics
#'
#' @keywords internal
register_built_in_metrics <- function() {
  register_metric(rmse_metric(), overwrite = TRUE)
  register_metric(bias_metric(), overwrite = TRUE)
  register_metric(coverage_metric(), overwrite = TRUE)
  register_metric(posterior_mean_metric(), overwrite = TRUE)
}
