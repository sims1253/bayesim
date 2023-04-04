#' Title
#'
#' @param family_list
#'
#' @return
#' @export
#'
#' @examples
get_prefit <- function(fit_conf, stan_pars, compile_dir = NULL) {
  family <- brms_family_lookup(
    fit_conf$fit_family,
    fit_conf$fit_link
  )
  if (is.null(fit_conf$prior)) {
    fit_conf$prior <- prior_lookup(fit_conf$fit_family)
  }
  formula <- brms::brmsformula(fit_conf$formula)
  data <- do.call(
    rng_lookup(fit_conf$fit_family),
    list(
      n = length(all.vars(formula$formula))
    )
  )
  names(data) <- all.vars(formula$formula)

  if(!is.null(compile_dir)) {
    prefit_name <- paste(
      fit_conf$fit_family, fit_conf$fit_link, fit_conf$formula, stan_pars$backend,
      sep = "_")
    prefit_dir <- paste0(
      paste(compile_dir, prefit_name, sep = "/"), ".RDS")
    prefit_dir <- gsub(" ", "", prefit_dir)

    if(file.exists(prefit_dir)) {
      cat("Found pre-compiled model", prefit_name, "update with new data\n")
      prefit <- readRDS(prefit_dir)
      prefit <- stats::update(prefit,
                              newdata = as.list(data),
                              formula. = formula,
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
  }

  # no pre-compiled data found, or feature unused, compile model now
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
  if(!is.null(compile_dir))
  {
    cat("Saved pre-compiled model", prefit_name, "\n")
    saveRDS(prefit, prefit_dir)
  }
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
build_prefit_list <- function(fit_configuration, stan_pars, compile_dir = NULL) {
  if (is.null(fit_configuration$prior)) {
    prefit_configurations <- unique(
      fit_configuration[c("fit_family", "fit_link", "formula")]
    )
  } else {
    prefit_configurations <- unique(
      fit_configuration[c("fit_family", "fit_link", "prior", "formula")]
    )
  }

  prefit_configurations <- lapply(split(
    prefit_configurations,
    sort(as.numeric(rownames(prefit_configurations)))
  ), as.list)

  prefit_list <- vector(mode = "list")
  for (conf in prefit_configurations) {
    prefit_list[[fit_conf_key(conf)]] <- get_prefit(
      fit_conf = conf, stan_pars, compile_dir
    )
  }

  return(prefit_list)
}
