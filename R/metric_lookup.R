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
                          vars_of_interest = NULL,
                          references = NULL,
                          threshold = 0.7,
                          psis_object = NULL,
                          ppred = NULL,
                          quantiles = NULL,
                          data_gen_output = NULL,
                          fit_conf = NULL,
                          ...) {
  # TODO this is a workaround. These shouldn't be lists in the first place.
  if (is.null(vars_of_interest)) {
    vars_of_interest <- data_gen_output$vars_of_interest
  } else {
    vars_of_interest <- unlist(vars_of_interest)
  }
  if (is.null(references)) {
    references <- data_gen_output$references
  } else {
    references <- unlist(references)
  }
  quantiles <- unlist(quantiles)
  y <- brms::get_y(fit)

  tryCatch(
    expr = {
      return(
        switch(metric,
          # Variable summaries
          "v_mean" = padd_variable_summay(
            draws, vars_of_interest, mean, metric
          ),
          "v_sd" = padd_variable_summay(
            draws, vars_of_interest, sd, metric
          ),
          "v_median" = padd_variable_summay(
            draws, vars_of_interest, median, metric
          ),
          "v_mad" = padd_variable_summay(
            draws, vars_of_interest, mad, metric
          ),
          "v_pos_prob" = padd_variable_summay(
            draws, vars_of_interest, bayeshear::variable_pos_prob, metric
          ),
          "v_quantiles" = padd_quantiles(draws, vars_of_interest, quantiles),

          # Variable distance measures
          "v_bias" = padd_variable_distance(
            draws, vars_of_interest, references, bayeshear::variable_bias, metric
          ),
          "v_rmse" = padd_variable_distance(
            draws, vars_of_interest, references, bayeshear::variable_rmse, metric
          ),
          "v_mae" = padd_variable_distance(
            draws, vars_of_interest, references, bayeshear::variable_mae, metric
          ),
          "v_mse" = padd_variable_distance(
            draws, vars_of_interest, references, bayeshear::variable_mse, metric
          ),
          "v_true_percentile" = padd_variable_distance(
            draws,
            vars_of_interest,
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
            list(
              pareto_k_values = list(bayeshear::pareto_k_values(psis_object))
            )
          },
          "time_per_sample" = bayeshear::sampling_time(fit, absolute = FALSE),

          # Variable MCMC Diagnostics
          "rhat" = {
            tmp <- as.list(posterior::rhat(fit, pars = vars_of_interest))
            names(tmp) <- lapply(vars_of_interest, function(x) paste0("rhat_", x))
            tmp
          },
          "ess_bulk" = {
            ess_list <- lapply(vars_of_interest, get_ess, fit, posterior::ess_bulk)
            names(ess_list) <- lapply(
              vars_of_interest, function(x) paste0("ess_bulk_", x)
            )
            ess_list
          },
          "ess_tail" = {
            ess_list <- lapply(vars_of_interest, get_ess, fit, posterior::ess_tail)
            names(ess_list) <- lapply(
              vars_of_interest, function(x) paste0("ess_tail_", x)
            )
            ess_list
          },

          # Predictive Metrics
          "elpd_loo" = elpd_loo_handler(fit),
          "elpd_loo_pointwise" = {
            loo_object <- brms::loo(fit)
            list(
              elpd_loo_pointwise = list(loo_object$pointwise[, 1]),
              mcse_elpd_loo_pointwise = list(loo_object$pointwise[, 2]),
              p_loo_pointwise = list(loo_object$pointwise[, 3])
            )
          },
          "elpd_loo_pointwise_summary" = elpd_pointwise_summaries(
            fit,
            quantiles
          ),
          "elpd_test" = elpd_test(fit, testing_data, FALSE),
          "elpd_test_pointwise_summary" =
            elpd_pointwise_summaries(fit, quantiles, testing_data),
          "rmse_loo" = rmse_loo(fit, psis_object = psis_object, yrep = ppred, ...),
          "rmse_loo_pointwise" = {
            loo_object <- rmse_loo(fit, psis_object = psis_object, yrep = ppred, return_object = TRUE)
            list(
              rmse_loo_pointwise = list(loo_object$pointwise[, 1]),
              rmse_loo_se = loo_object$estimates[1, 2]
            )
          },
          "rmse_loo_pointwise_summary" =
            get_custom_loo_summary(
              rmse_loo(fit,
                psis_object,
                yrep = ppred,
                return_object = TRUE
              ),
              quantiles, "rmse_loo"
            ),
          "rmse_test" = rmse_test(fit, testing_data),
          "rmse_test_pointwise_summary" =
            get_custom_loo_summary(
              rmse_test(fit,
                testing_data,
                return_object = TRUE
              ),
              quantiles, "rmse_test"
            ),
          "r2_loo" = r2_loo(fit, psis_object = psis_object, yrep = ppred),
          "r2_loo_pointwise" = {
            loo_object <- r2_loo(fit, psis_object = psis_object, yrep = ppred, return_object = TRUE)
            list(
              r2_loo_pointwise = list(loo_object$pointwise[, 1]),
              r2_loo_se = loo_object$estimates[1, 2]
            )
          },
          "r2_loo_pointwise_summary" =
            get_custom_loo_summary(
              r2_loo(fit,
                psis_object,
                return_object = TRUE
              ),
              quantiles, "r2_loo"
            ),
          "r2_test" = r2_test(fit, testing_data),
          "r2_test_pointwise_summary" =
            get_custom_loo_summary(
              r2_test(fit,
                testing_data,
                return_object = TRUE
              ),
              quantiles, "r2_test"
            ),

          # Posterior sample based metrics
          "log_lik_pointwise" = {
            ll <- brms::log_lik(fit)
            list(
              log_lik_pointwise_mean = list(colMeans(ll)),
              log_lik_pointwise_sd = list(apply(ll, 2, sd))
            )
          },
          "log_lik_summary" = observation_x_sample_summarizer(
            brms::log_lik(fit),
            quantiles,
            "log_lik_summary"
          ),
          "ppred_summary_y_scaled" = observation_x_sample_summarizer(
            ((ppred - mean(y)) / sd(y)),
            quantiles,
            "ppred_summary"
          ),
          "ppred_pointwise" = {
            list(
              ppred_pointwise_mean = list(colMeans(ppred)),
              ppred_pointwise_sd = list(apply(ppred, 2, sd))
            )
          },
          "residuals" =
            list(residuals = list(residuals(fit, method = "posterior_predict")[, 1])),
          "posterior_linpred" = {
            linpred <- brms::posterior_linpred(fit)
            list(
              posterior_linpred_mean = list(colMeans(linpred)),
              posterior_linpred_sd = list(apply(linpred, 2, sd))
            )
          },
          "posterior_linpred_transformed" = {
            linpred <- do.call(
              link_lookup(fit_conf$fit_link, inv = TRUE),
              list(brms::posterior_linpred(fit))
            )
            list(
              posterior_linpred_transformed_mean = list(colMeans(linpred)),
              posterior_linpred_transformed_sd = list(apply(linpred, 2, sd))
            )
          },
          # Observations
          "y_pointwise" = {
            list(y_pointwise = list(y))
          },
          "y_pointwise_z_scaled" = {
            tmp <- (as.list(y) - mean(y)) / sd(y)
            names(tmp) <- lapply(
              seq_along(tmp),
              function(x) paste0("obs_z_scaled", x)
            )
            tmp
          },
          "y_summaries" = list(
            y_mean = mean(y),
            y_sd = sd(y)
          ),

          # Data
          "data_gen" = data_gen_output,

          # Fits
          "fit_gen" = fit_conf,
          stop(paste(metric, "is not a supported metric!"))
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
  quantile_list[paste0(name, "_mean")] <- mean(pointwise)
  quantile_list[paste0(name, "_sd")] <- sd(pointwise)
  quantile_list[paste0(name, "_se_mean")] <- loo_object$estimates[1, 2] /
    length(pointwise)
  return(quantile_list)
}
