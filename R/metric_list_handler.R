#' Collect multiple metrics at once
#'
#' A convenience function that collects all given metrics at once. This
#' can save time compared to manually calling all metric functions
#' individually, as some variables can be reused instead of being calculated
#' multiple times.
#'
#' @param fit A brmsfit object.
#' @param metrics A vector of metric identifiers. See details for supported
#' identifiers.
#'
#' @return A named list containing all requested metrics' results.
#'
#' @param data_gen_output Output from data generation containing true parameters
#' @param fit_conf Model configuration list or data.frame row
#' @param ... Additional arguments passed to metric calculation functions
#'
#' @details Currently, the following identifiers are supported:
#' \itemize{
#' \item "bias": posterior bias
#' \item "divergents": number of divergent transitions
#' \item "ess": effective sample size
#' \item "elpd_loo": ELPD-LOO
#' \item "elpd_newdata": ELPD on new data
#' \item "epred": expected predictions
#' \item "mae_s": mean absolute error scaled
#' \item "p_mean": posterior mean
#' \item "p_sd": posterior standard deviation
#' \item "pareto_k": Pareto k diagnostic values
#' \item "pos_prob": posterior probability of positive values
#' \item "ppred": posterior predictions
#' \item "pq": posterior quantiles
#' \item "q_true": true quantiles
#' \item "r2_loo": LOO R-squared
#' \item "r2_newdata": R-squared on new data
#' \item "residuals": model residuals
#' \item "rhat": R-hat diagnostic
#' \item "rmse_loo": LOO root mean square error
#' \item "rmse_newdata": RMSE on new data
#' \item "rmse_s": scaled RMSE
#' \item "rstar": R* diagnostic
#' \item "time_sampling": sampling time
#' \item "time_total": total time
#' \item "time_warmup": warmup time
#' \item "y": observed values
#' }
#'
#' Note,that not all identifiers are supported for each input class.
#'
#' @export
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Requires brms package and a fitted model
#' fit <- brms::brm(y ~ x, data = mydata)
#' metric_list_handler(fit, c("elpd_loo", "rhat"), data_gen_output, fit_conf)
#' }
metric_list_handler <- function(fit, metrics, data_gen_output, fit_conf, ...) {
  needs_draws <- list(
    "v_mean",
    "v_sd",
    "v_median",
    "v_mad",
    "v_pos_prob",
    "quantiles",
    "v_bias",
    "v_rmse",
    "v_mae",
    "v_mse",
    "v_percentile",
    "rstar"
  )
  needs_psis <- list(
    "pareto_k_values",
    "bad_pareto_ks",
    "rmse_loo",
    "r2_loo",
    "rmse_loo_pointwise_summary",
    "r2_loo_pointwise",
    "r2_loo_pointwise_summary"
  )
  needs_ppred <- list(
    "rmse_loo",
    "r2_loo",
    "rmse_loo_pointwise_summary",
    "r2_loo_pointwise",
    "r2_loo_pointwise_summary",
    "ppred_pointwise",
    "ppred_summary_y_scaled"
  )
  precompute_error <- tryCatch(
    {
      if (length(intersect(needs_draws, metrics)) > 0) {
        draws <- posterior::as_draws(fit)
      } else {
        draws <- NULL
      }
      if (length(intersect(needs_psis, metrics)) > 0) {
        psis_object <- tryCatch(
          get(".psis", envir = asNamespace("brms"))(
            fit,
            newdata = fit$data,
            resp = NULL
          ),
          error = function(e) NULL
        )
      } else {
        psis_object <- NULL
      }
      if (length(intersect(needs_ppred, metrics)) > 0) {
        ppred <- brms::posterior_predict(fit, fit$data)
      } else {
        ppred <- NULL
      }

      NULL
    },
    error = function(e) {
      dplyr::as_tibble(do.call(c, list(data_gen_output, fit_conf)))
    }
  )

  if (!is.null(precompute_error)) {
    return(precompute_error)
  }

  results <- lapply(
    metrics,
    metric_lookup,
    fit = fit,
    draws = draws,
    psis_object = psis_object,
    ppred = ppred,
    data_gen_output = data_gen_output,
    fit_conf = fit_conf,
    ...
  )
  dplyr::as_tibble(do.call(c, results))
}
