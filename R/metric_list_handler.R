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
#' @details Currently, the following identifiers are supported. See linked
#'  functions for required additional arguments:
#' \itemize{
#' \item "bias": [posterior_bias()]
#' \item "divergents": [divergents()]
#' \item "ess": [ess()]
#' \item "elpd_loo": [elpd_loo()]
#' \item "elpd_newdata": [elpd_newdata()]
#' \item "epred": [epred()]
#' \item "mae_s": \code{\link{}}
#' \item "p_mean": [p_mean()]
#' \item "p_sd": [p_sd()]
#' \item "pareto_k": \code{\link{}}
#' \item "pos_prob": \code{\link{}}
#' \item "ppred": \code{\link{}}
#' \item "pq": \code{\link{}}
#' \item "q_true": \code{\link{}}
#' \item "r2_loo": \code{\link{}}
#' \item "r2_newdata": \code{\link{}}
#' \item "residuals": \code{\link{}}
#' \item "rhat": \code{\link{}}
#' \item "rmse_loo": \code{\link{}}
#' \item "rmse_newdata": \code{\link{}}
#' \item "rmse_s": \code{\link{}}
#' \item "rstar": \code{\link{}}
#' \item "time_sampling": \code{\link{}}
#' \item "time_total": \code{\link{}}
#' \item "time_warmup": \code{\link{}}
#' \item "y": \code{\link{}}
#' }
#'
#' Note,that not all identifiers are supported for each input class.
#'
#' @export
#' @examples
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
  result <- tryCatch(
    {
      if (length(intersect(needs_draws, metrics)) > 0) {
        draws <- posterior::as_draws(fit)
      } else {
        draws <- NULL
      }
      if (length(intersect(needs_psis, metrics)) > 0) {
        psis_object <- brms:::.psis(fit, newdata = fit$data, resp = NULL)
      } else {
        psis_object <- NULL
      }
      if (length(intersect(needs_ppred, metrics)) > 0) {
        ppred <- brms::posterior_predict(fit, fit$data)
      } else {
        ppred <- NULL
      }
    },
    error = function(e) {
      return(dplyr::as_tibble(do.call(c, list(data_gen_output, fit_conf))))
    }
  )

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
  return(dplyr::as_tibble(do.call(c, results)))
}
