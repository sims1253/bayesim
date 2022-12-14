#' Access metrics with string identifiers
#'
#' This function is mostly a helper function that maps identifier strings to
#' metrics for convenient use in some contexts.
#'
#' @param metric A string that identifies a supported metric
#'
#' @return The function corresponding to the identifier string.
#' @export
#'
metric_lookup <- function(metric,
                          fit = NULL,
                          draws = NULL,
                          testing_data = NULL,
                          variables = NULL,
                          references = NULL,
                          threshold = 0.7,
                          psis_object = NULL,
                          quantiles = NULL,
                          ...) {
  # TODO this is a workaround. These shouldn't be lists in the first place.
  variables <- unlist(variables)
  references <- unlist(references)
  quantiles <- unlist(quantiles)

  tryCatch(
    expr = {
      return(
        switch(metric,
          # Variable summaries
          "v_mean" = padd_variable_summay(
            draws, variables, mean, metric
          ),
          "v_sd" = padd_variable_summay(
            draws, variables, sd, metric
          ),
          "v_median" = padd_variable_summay(
            draws, variables, median, metric
          ),
          "v_mad" = padd_variable_summay(
            draws, variables, mad, metric
          ),
          "v_pos_prob" = padd_variable_summay(
            draws, variables, bayeshear::variable_pos_prob, metric
          ),
          "v_quantiles" = padd_quantiles(draws, variables, quantiles),

          # Variable distance measures
          "v_bias" = padd_variable_distance(
            draws, variables, references, bayeshear::variable_bias, metric
          ),
          "v_rmse" = padd_variable_distance(
            draws, variables, references, bayeshear::variable_rmse, metric
          ),
          "v_mae" = padd_variable_distance(
            draws, variables, references, bayeshear::variable_mae, metric
          ),
          "v_mse" = padd_variable_distance(
            draws, variables, references, bayeshear::variable_mse, metric
          ),
          "v_true_percentile" = padd_variable_distance(
            draws,
            variables,
            references,
            bayeshear::variable_true_percentile,
            metric
          ),

          # Global MCMC Diagnostics
          "divergent_transitions_rel" = list(
            "divergent_transitions_rel" = bayeshear::divergent_transitions(fit)
          ),
          "divergent_transitions_abs" =
            list(
              "divergent_transitions_abs" =
                bayeshear::divergent_transitions(fit, absolute = TRUE)
            ),
          "rstar" = list("rstar" = posterior::rstar(draws)),
          "bad_pareto_ks" = list(
            "bad_pareto_ks" =
              bayeshear::bad_pareto_ks(fit, threshold)
          ),
          "pareto_k_values" = {
            tmp <- as.list(bayeshear::pareto_k_values(psis_object))
            names(tmp) <- lapply(
              seq_along(tmp), function(x) paste0("pareto_k_obs_", x))
            tmp
          },
          "time" = bayeshear::sampling_time(fit, absolute = FALSE),

          # Variable MCMC Diagnostics
          "rhat" = {
            tmp <- as.list(posterior::rhat(fit, pars = variables))
            names(tmp) <- lapply(variables, function(x) paste0("rhat_", x))
            tmp
          },
          "ess_bulk" = {
            ess_list <- lapply(variables, get_ess, fit, posterior::ess_bulk)
            names(ess_list) <- lapply(
              variables, function(x) paste0("ess_bulk_", x)
            )
            ess_list
          },
          "ess_tail" = {
            ess_list <- lapply(variables, get_ess, fit, posterior::ess_tail)
            names(ess_list) <- lapply(
              variables, function(x) paste0("ess_tail_", x)
            )
            ess_list
          },

          # Predictive Metrics
          "elpd_loo_pointwise_summary" = elpd_pointwise_summaries(
            fit,
            quantiles
          ),
          "elpd_loo" = elpd_loo_handler(fit),
          "elpd_test" = elpd_test(fit, testing_data, FALSE),
          "elpd_test_pointwise_summary" =
            elpd_pointwise_summaries(fit, quantiles, testing_data),
          "rmse_loo" = rmse_loo(fit, ...),
          "rmse_loo_pointwise_summary" =
            get_custom_loo_summary(rmse_loo(fit,
                                            psis_object,
                                            return_object = TRUE),
                                   quantiles, "rmse_loo"),
          "rmse_test" = rmse_test(fit, testing_data),
          "rmse_test_pointwise_summary" =
            get_custom_loo_summary(rmse_test(fit,
                                             testing_data,
                                            return_object = TRUE),
                                   quantiles, "rmse_test"),
          "r2_loo" = r2_loo(fit, ...),
          "r2_loo_pointwise_summary" =
            get_custom_loo_summary(r2_loo(fit,
                                          psis_object,
                                          return_object = TRUE),
                                   quantiles, "r2_loo"),
          "r2_test" = r2_test(fit, testing_data),
          "r2_test_pointwise_summary" =
            get_custom_loo_summary(r2_test(fit,
                                           testing_data,
                                           return_object = TRUE),
                                   quantiles, "r2_test"),

          # Posterior sample based metrics
          "log_lik_summary" = observation_x_sample_summarizer(
            brms::log_lik(fit),
            quantiles,
            "log_lik_summary"
          ),
          "ppred_summary" = observation_x_sample_summarizer(
            brms::posterior_predict(fit),
            quantiles,
            "ppred_summary"
          ),
          "residuals" =
            list(residuals = residuals(fit, method = "posterior_predict")[, 1]),

          # Observations
          "y" = {
            tmp <- as.list(brms::get_y(fit))
            names(tmp) <- lapply(seq_along(tmp), function(x) paste0("obs_", x))
            tmp
          },
          "y_summaries" = list(
            y_mean = mean(brms::get_y(fit)),
            y_sd = sd(brms::get_y(fit))
          )
        )
      )
    }, error = function(e) {
      return(list())
    }
  )
}

padd_variable_summay <- function(draws, variables, metric, name) {
  tmp <- as.list(bayeshear::variable_summary(draws, variables, metric))
  names(tmp) <- lapply(
    variables,
    function(x) paste0(name, "_", x)
  )
  return(tmp)
}

padd_variable_distance <- function(draws, variables, references, metric, name) {
  tmp <- as.list(
    bayeshear::variable_distance(draws, variables, references, metric)
  )
  names(tmp) <- lapply(
    variables,
    function(x) paste0(name, "_", x)
  )
  return(tmp)
}

padd_quantiles <- function(draws, variables, quantiles) {
  quantile_list <- bayeshear::posterior_quantiles(draws, variables, quantiles)
  quantile_list <- unlist(quantile_list, recursive = FALSE)
  names(quantile_list) <- gsub("[.]", "_", names(quantile_list))
  names(quantile_list) <- gsub("[%]", "pq", names(quantile_list))
  return(as.list(quantile_list))
}

get_ess <- function(variable, fit, fun) {
  return(
    do.call(
      fun,
      list(posterior::extract_variable_matrix(fit, variable = variable))
    )
  )
}

observation_x_sample_summarizer <- function(sample_matrix, quantiles, name) {
  out <- vector(mode = "list", length = ncol(sample_matrix))

  for (i in seq_len(ncol(sample_matrix))) {
    tmp <- as.list(quantile(sample_matrix[, i], prob = quantiles))
    names(tmp) <- lapply(
      quantiles,
      function(x) paste0(name, "_obs_", i, "_quantile_", x)
    )
    out[[i]] <- tmp
    out[[i]][paste0(name, "_obs_", i, "_mean")] <- mean(sample_matrix[, i])
    out[[i]][paste0(name, "_obs_", i, "_sd")] <- sd(sample_matrix[, i])
  }
  return(unlist(out, recursive = FALSE))
}

get_col_summaries <- function(column, quantiles, name) {
  tmp <- as.list(quantile(column, prob = quantiles))
  names(tmp) <- lapply(
    quantiles,
    function(x) paste0(name, "_obs_", i, "_quantile_", x)
  )
  out[[i]] <- tmp
  out[[i]][paste0(name, "_obs_", i, "_mean")] <- mean(column)
  out[[i]][paste0(name, "_obs_", i, "_sd")] <- sd(column)
  return(out)
}

get_custom_loo_summary <- function(loo_object, quantiles, name) {
  pointwise <- loo_object$pointwise
  quantile_list <- as.list(quantile(pointwise, prob = quantiles))
  names(quantile_list) <- lapply(
    quantiles,
    function(x) paste0(name, "_quantile_", x)
  )
  out[paste0(name, "_mean")] <- mean(pointwise)
  out[paste0(name, "_sd")] <- sd(pointwise)
  out[paste0(name, "_se_mean")] <- loo_object$estimates[1, 2] /
                                   length(pointwise)
  return(out)
}
