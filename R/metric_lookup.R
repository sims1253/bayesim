#' Title
#'
#' @param identifier
#' @param fit
#' @param ...
#' @param posterior_draws
#'
#' @return
#' @export
#'
#' @examples
metric_lookup <- function(identifier, fit, posterior_draws = NULL, ...) {
  if (is.null(posterior_draws)) {
    posterior_draws <- as.vector(posterior::as_draws_array(fit, variable = "b_x"))
  }

  if (grepl("^hyp_", identifier)) {
    hyp <- substring(identifier, 5)
    return(hypothesis_cis(hyp, fit, ...))
  }
  if (grepl("pq_", identifier)) {
    prob <- as.numeric(substring(identifier, 4))
    return(p_quantiles(posterior_draws, prob))
  } else {
    switch(identifier,
      "divergents" = return(divergents(fit)),
      "rstar" = return(p_rstar(fit)),
      "pareto_k" = return(bad_pareto_ks(fit)),
      "time" = return(sampling_time(fit)),
      "rhat" = return(brms::rhat(fit)["b_x"]),
      "ess_bulk" = return(rstan::ess_bulk(posterior_draws)),
      "ess_tail" = return(rstan::ess_tail(posterior_draws)),
      "q_true" = return(q_true(posterior_draws, ...)),
      "bias" = return(p_bias(posterior_draws, ...)),
      "rmse_s" = return(rmse_s(posterior_draws, ...)),
      "mae_s" = return(mae_s(posterior_draws, ...)),
      "p_mean" = return(p_mean(posterior_draws)),
      "p_sd" = return(p_sd(posterior_draws))
    )
  }
}
