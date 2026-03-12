#' Access metrics with string identifiers
#'
#' This function is mostly a helper function that maps identifier strings to
#' metrics for convenient use in some contexts.
#'
#' @param metric A string that identifies a supported metric
#' @param fit A brms fit object
#' @param draws Posterior draws object
#' @param testing_data Data used for testing/out-of-sample evaluation
#' @param vars_of_interest Variables to compute metrics for
#' @param references Reference values for computing distance metrics
#' @param threshold Threshold for diagnostic checks (default: 0.7)
#' @param psis_object PSIS-LOO object for importance sampling
#' @param ppred Posterior predictive samples matrix
#' @param quantiles Quantiles to compute for summaries
#' @param data_gen_output Output from data generation process
#' @param fit_conf Fit configuration list
#' @param ... Additional arguments passed to metric functions
#'
#' @return The function corresponding to the identifier string.
#' @export
#'
metric_lookup <- function(
  metric,
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
  ...
) {
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
        switch(
          metric,
          # Variable summaries
          "v_mean" = padd_variable_summay(
            draws,
            vars_of_interest,
            mean,
            metric
          ),
          "v_sd" = padd_variable_summay(
            draws,
            vars_of_interest,
            sd,
            metric
          ),
          "v_median" = padd_variable_summay(
            draws,
            vars_of_interest,
            median,
            metric
          ),
          "v_mad" = padd_variable_summay(
            draws,
            vars_of_interest,
            mad,
            metric
          ),
          "v_pos_prob" = padd_variable_summay(
            draws,
            vars_of_interest,
            function(x) {
              tryCatch(bayeshear::variable_pos_prob(x), error = function(e) {
                NA_real_
              })
            },
            metric
          ),
          "v_quantiles" = padd_quantiles(draws, vars_of_interest, quantiles),

          # Variable distance measures
          "v_bias" = padd_variable_distance(
            draws,
            vars_of_interest,
            references,
            function(d, v, r) {
              tryCatch(bayeshear::variable_bias(d, v, r), error = function(e) {
                NA_real_
              })
            },
            metric
          ),
          "v_rmse" = padd_variable_distance(
            draws,
            vars_of_interest,
            references,
            function(d, v, r) {
              tryCatch(bayeshear::variable_rmse(d, v, r), error = function(e) {
                NA_real_
              })
            },
            metric
          ),
          "v_mae" = padd_variable_distance(
            draws,
            vars_of_interest,
            references,
            function(d, v, r) {
              tryCatch(bayeshear::variable_mae(d, v, r), error = function(e) {
                NA_real_
              })
            },
            metric
          ),
          "v_mse" = padd_variable_distance(
            draws,
            vars_of_interest,
            references,
            function(d, v, r) {
              tryCatch(bayeshear::variable_mse(d, v, r), error = function(e) {
                NA_real_
              })
            },
            metric
          ),
          "v_true_percentile" = padd_variable_distance(
            draws,
            vars_of_interest,
            references,
            function(d, v, r) {
              tryCatch(
                bayeshear::variable_true_percentile(d, v, r),
                error = function(e) NA_real_
              )
            },
            metric
          ),

          # Global MCMC Diagnostics
          "divergent_transitions_rel" = list(
            "divergent_transitions_rel" = tryCatch(
              bayeshear::divergent_transitions(fit),
              error = function(e) NA_real_
            )
          ),
          "divergent_transitions_abs" = list(
            "divergent_transitions_abs" = tryCatch(
              bayeshear::divergent_transitions(fit, absolute = TRUE),
              error = function(e) NA_real_
            )
          ),
          "rstar" = list(
            "rstar" = tryCatch(posterior::rstar(draws), error = function(e) {
              NA_real_
            })
          ),
          "bad_pareto_ks" = list(
            "bad_pareto_ks" = tryCatch(
              bayeshear::bad_pareto_ks(fit, threshold),
              error = function(e) NA_real_
            )
          ),
          "pareto_k_values" = {
            list(
              pareto_k_values = list(tryCatch(
                bayeshear::pareto_k_values(psis_object),
                error = function(e) list()
              ))
            )
          },
          "time_per_sample" = tryCatch(
            bayeshear::sampling_time(fit, absolute = FALSE),
            error = function(e) NA_real_
          ),

          # Variable MCMC Diagnostics
          "rhat" = {
            tmp <- as.list(posterior::rhat(fit, pars = vars_of_interest))
            names(tmp) <- lapply(vars_of_interest, function(x) {
              paste0("rhat_", x)
            })
            tmp
          },
          "ess_bulk" = {
            ess_list <- lapply(
              vars_of_interest,
              get_ess,
              fit,
              posterior::ess_bulk
            )
            names(ess_list) <- lapply(
              vars_of_interest,
              function(x) paste0("ess_bulk_", x)
            )
            ess_list
          },
          "ess_tail" = {
            ess_list <- lapply(
              vars_of_interest,
              get_ess,
              fit,
              posterior::ess_tail
            )
            names(ess_list) <- lapply(
              vars_of_interest,
              function(x) paste0("ess_tail_", x)
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
          "elpd_test_pointwise_summary" = elpd_pointwise_summaries(
            fit,
            quantiles,
            testing_data
          ),
          "rmse_loo" = rmse_loo(
            fit,
            psis_object = psis_object,
            yrep = ppred,
            ...
          ),
          "rmse_loo_pointwise" = {
            loo_object <- rmse_loo(
              fit,
              psis_object = psis_object,
              yrep = ppred,
              return_object = TRUE
            )
            list(
              rmse_loo_pointwise = list(loo_object$pointwise[, 1]),
              rmse_loo_se = loo_object$estimates[1, 2]
            )
          },
          "rmse_loo_pointwise_summary" = get_custom_loo_summary(
            rmse_loo(fit, psis_object, yrep = ppred, return_object = TRUE),
            quantiles,
            "rmse_loo"
          ),
          "rmse_test" = rmse_test(fit, testing_data),
          "rmse_test_pointwise_summary" = get_custom_loo_summary(
            rmse_test(fit, testing_data, return_object = TRUE),
            quantiles,
            "rmse_test"
          ),
          "r2_loo" = r2_loo(fit, psis_object = psis_object, yrep = ppred),
          "r2_loo_pointwise" = {
            loo_object <- r2_loo(
              fit,
              psis_object = psis_object,
              yrep = ppred,
              return_object = TRUE
            )
            list(
              r2_loo_pointwise = list(loo_object$pointwise[, 1]),
              r2_loo_se = loo_object$estimates[1, 2]
            )
          },
          "r2_loo_pointwise_summary" = get_custom_loo_summary(
            r2_loo(fit, psis_object, return_object = TRUE),
            quantiles,
            "r2_loo"
          ),
          "r2_test" = r2_test(fit, testing_data),
          "r2_test_pointwise_summary" = get_custom_loo_summary(
            r2_test(fit, testing_data, return_object = TRUE),
            quantiles,
            "r2_test"
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
          "residuals" = list(
            residuals = list(residuals(fit, method = "posterior_predict")[, 1])
          ),
          "posterior_linpred" = {
            linpred <- brms::posterior_linpred(fit)
            list(
              posterior_linpred_mean = list(colMeans(linpred)),
              posterior_linpred_sd = list(apply(linpred, 2, sd))
            )
          },
          "posterior_linpred_transformed" = {
            linpred <- tryCatch(
              do.call(
                bayesfam::link_lookup(fit_conf$fit_link, inv = TRUE),
                list(brms::posterior_linpred(fit))
              ),
              error = function(e) brms::posterior_linpred(fit)
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
    },
    error = function(e) {
      list()
    }
  )
}

padd_variable_summay <- function(draws, variables, metric, name) {
  tmp <- tryCatch(
    as.list(bayeshear::variable_summary(draws, variables, metric)),
    error = function(e) {
      # Fallback if bayeshear is not available
      tmp <- lapply(variables, function(v) {
        vals <- posterior::extract_variable_matrix(draws, variable = v)
        metric(vals)
      })
      names(tmp) <- variables
      tmp
    }
  )
  names(tmp) <- lapply(
    variables,
    function(x) paste0(name, "_", x)
  )
  tmp
}

padd_variable_distance <- function(draws, variables, references, metric, name) {
  tmp <- tryCatch(
    as.list(
      bayeshear::variable_distance(draws, variables, references, metric)
    ),
    error = function(e) {
      # Return NA for all variables if bayeshear is not available
      tmp <- rep(NA_real_, length(variables))
      names(tmp) <- variables
      as.list(tmp)
    }
  )
  names(tmp) <- lapply(
    variables,
    function(x) paste0(name, "_", x)
  )
  tmp
}

padd_quantiles <- function(draws, variables, quantiles) {
  quantile_list <- tryCatch(
    bayeshear::posterior_quantiles(draws, variables, quantiles),
    error = function(e) {
      # Fallback using posterior package
      lapply(variables, function(v) {
        vals <- posterior::extract_variable_matrix(draws, variable = v)
        quantile(vals, probs = quantiles)
      })
    }
  )
  quantile_list <- unlist(quantile_list, recursive = FALSE)
  names(quantile_list) <- gsub("[.]", "_", names(quantile_list))
  names(quantile_list) <- gsub("[%]", "pq", names(quantile_list))
  as.list(quantile_list)
}

get_ess <- function(variable, fit, fun) {
  do.call(
    fun,
    list(posterior::extract_variable_matrix(fit, variable = variable))
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
  unlist(out, recursive = FALSE)
}

#' Compute column summaries with quantiles
#'
#' @param column A numeric vector to summarize
#' @param quantiles Numeric vector of quantile probabilities
#' @param name A string prefix for output names
#'
#' @return A named list with quantile values, mean, and standard deviation
#' @keywords internal
get_col_summaries <- function(column, quantiles, name) {
  tmp <- as.list(quantile(column, prob = quantiles))
  names(tmp) <- lapply(
    quantiles,
    function(x) paste0(name, "_quantile_", x)
  )
  tmp[[paste0(name, "_mean")]] <- mean(column)
  tmp[[paste0(name, "_sd")]] <- sd(column)
  tmp
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
  quantile_list
}
