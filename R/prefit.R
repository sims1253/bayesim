#' Create a precompiled brms model object
#'
#' Creates a brmsfit object with compiled Stan code that can be efficiently
#' updated with new data during simulation runs, avoiding recompilation.
#'
#' @param fit_conf A list or data.frame row containing model configuration:
#'   \itemize{
#'     \item fit_family: The model family (e.g., "gaussian", "student_t")
#'     \item fit_link: The link function (e.g., "identity", "log")
#'     \item formula: The model formula as a character string
#'     \item prior: (optional) Prior specifications
#'   }
#' @param stan_pars A list containing Stan parameters, must include:
#'   \itemize{
#'     \item backend: Stan backend ("cmdstanr" or "rstan")
#'   }
#'
#' @return A brmsfit object with compiled Stan code (chains = 0)
#' @export
#'
#' @examples
#' \dontrun{
#' fit_conf <- list(
#'   fit_family = "gaussian",
#'   fit_link = "identity",
#'   formula = "y ~ x"
#' )
#' stan_pars <- list(backend = "cmdstanr")
#' prefit <- get_prefit(fit_conf, stan_pars)
#' }
get_prefit <- function(fit_conf, stan_pars) {
  family <- tryCatch(
    bayesfam::brms_family_lookup(fit_conf$fit_family, fit_conf$fit_link),
    error = function(e) {
      # Fallback to brms family lookup if bayesfam is not available
      brms::brmsfamily(fit_conf$fit_family, fit_conf$fit_link)
    }
  )
  if (is.null(fit_conf$prior)) {
    fit_conf$prior <- prior_lookup(fit_conf$fit_family)
  }
  formula <- brms::brmsformula(fit_conf$formula)
  data <- tryCatch(
    do.call(
      bayesfam::rng_lookup(fit_conf$fit_family),
      list(n = length(all.vars(formula$formula)))
    ),
    error = function(e) {
      # Fallback: generate simple data
      n_vars <- length(all.vars(formula$formula))
      set.seed(12345)
      as.data.frame(matrix(rnorm(n_vars * 10), ncol = n_vars))
    }
  )
  names(data) <- all.vars(formula$formula)
  prefit <- brms::brm(
    formula = formula,
    data = as.list(data),
    family = family,
    stanvars = family$stanvars,
    chains = 0,
    refresh = 0,
    silent = 2,
    backend = stan_pars$backend,
    prior = fit_conf$prior,
    init = 0.1
  )
  prefit
}


#' Prepare a list of brmsfit objects to update during repeated simulations
#'
#' This is mainly used to save on constant recompilation times.
#'
#' @param fit_configuration A named list that currently holds family, link and
#'                          prior.
#' @param stan_pars A named list which contains a backend field.
#'
#' @return A named list of precompiled fit objects keyed by fit configuration.
#' @export
build_prefit_list <- function(fit_configuration, stan_pars) {
  if (is.null(fit_configuration$prior)) {
    prefit_configurations <- unique(
      fit_configuration[c("fit_family", "fit_link", "formula")]
    )
  } else {
    prefit_configurations <- unique(
      fit_configuration[c("fit_family", "fit_link", "prior", "formula")]
    )
  }

  prefit_configurations <- lapply(
    split(
      prefit_configurations,
      sort(as.numeric(rownames(prefit_configurations)))
    ),
    as.list
  )

  prefit_list <- vector(mode = "list")
  for (conf in prefit_configurations) {
    prefit_list[[fit_conf_key(conf)]] <- get_prefit(
      fit_conf = conf,
      stan_pars
    )
  }

  prefit_list
}
