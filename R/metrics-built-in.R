#' @title RMSE Metric
#' @description Root Mean Square Error between predictions and true values.
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
  register_metric(RmseMetric(name = "rmse"), overwrite = TRUE)
  register_metric(BiasMetric(name = "bias"), overwrite = TRUE)
  register_metric(CoverageMetric(name = "coverage"), overwrite = TRUE)
  register_metric(
    PosteriorMeanMetric(name = "posterior_mean"),
    overwrite = TRUE
  )
}
