#' Title
#'
#' @param family_list
#'
#' @return
#' @export
#'
#' @examples
get_prefit <- function(fit_conf, stan_pars) {
  family <- brms_family_lookup(
    fit_conf$fit_family,
    fit_conf$fit_link
  )
  if (is.null(fit_conf$prior)) {
    fit_conf$prior <- prior_lookup(fit_conf$fit_family)
  }
  prefit <- brms::brm(
    y ~ 1 + x,
    data = list(y = c(0.5), x = c(1)),
    family = family,
    stanvars = family$stanvars,
    chains = 0,
    refresh = 0,
    silent = 2,
    backend = stan_pars$backend,
    prior = fit_conf$prior,
    init = 0.1
  )
  return(prefit)
}


#' Prepare a list of brmsfit objects to update during repeated simulations
#'
#' This is mainly used to save on constant recompilation times.
#'
#' @param fit_configuration A named list that currently holds family, link and
#'                          prior.
#' @param stan_pars A named list which contains a backend field.
#'
#' @return A list of
#' @export
#'
#' @examples
build_prefit_list <- function(fit_configuration, stan_pars) {
  if (is.null(fit_configuration$prior)) {
    prefit_configurations <- unique(
      fit_configuration[c("fit_family", "fit_link")]
    )
  } else {
    prefit_configurations <- unique(
      fit_configuration[c("fit_family", "fit_link", "prior")]
    )
  }

  prefit_configurations <- lapply(split(
    prefit_configurations,
    sort(as.numeric(rownames(prefit_configurations)))
  ), as.list)

  prefit_list <- vector(mode = "list")
  for (conf in prefit_configurations) {
    prefit_list[[fit_conf_key(conf)]] <- get_prefit(
      fit_conf = conf, stan_pars
    )
  }

  return(prefit_list)
}
