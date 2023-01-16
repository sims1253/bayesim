#' Title
#'
#' @param family
#' @param link
#'
#' @return
#' @export
#'
#' @examples
brms_family_lookup <- function(family, link = NULL) {
  switch(family,
    "beta" = brms::brmsfamily("beta", link = link),
    "kumaraswamy" = bayesfam::kumaraswamy(link = link),
    "logitnormal" = bayesfam::logitnormal(link = link),
    "cauchitnormal" = bayesfam::cauchitnormal(link = link),
    "cloglognormal" = bayesfam::cloglognormal(link = link),
    "simplex" = bayesfam::simplex(link = link),
    "gaussian" = brms::brmsfamily("gaussian", link = link),
    "gamma" = brms::brmsfamily("gamma", link = link),
    "weibull" = brms::brmsfamily("weibull", link = link),
    "lognormal" = brms::brmsfamily("lognormal", link = link),
    "softplusnormal" = bayesfam::softplusnormal(link = link),
    "lomax" = bayesfam::lomax(link = link),
    "frechet" = brms::brmsfamily("frechet", link = link),
    "inverse.gaussian" = brms::brmsfamily("inverse.gaussian", link = link),
    "betaprime" = bayesfam::betaprime(link = link),
    "gompertz" = bayesfam::gompertz(link = link)
  )
}

#' Title
#'
#' @param family
#' @param link
#'
#' @return
#' @export
#'
#' @examples
rng_lookup <- function(family) {
  switch(family,
    "beta" = bayesfam::rbeta_mean,
    "kumaraswamy" = bayesfam::rkumaraswamy,
    "logitnormal" = bayesfam::rlogitnormal,
    "cauchitnormal" = bayesfam::rcauchitnormal,
    "cloglognormal" = bayesfam::rcloglognormal,
    "simplex" = bayesfam::rsimplex,
    "gaussian" = rnorm,
    "gamma" = bayesfam::rgamma_mean,
    "weibull" = bayesfam::rweibull_median,
    "lognormal" = bayesfam::rlognormal,
    "softplusnormal" = bayesfam::rsoftplusnormal,
    "lomax" = bayesfam::rlomax,
    "frechet" = bayesfam::rfrechet_median,
    "inverse.gaussian" = brms::rinv_gaussian,
    "betaprime" = bayesfam::rbetaprime,
    "gompertz" = bayesfam::rgompertz
  )
}


#' Title
#'
#' @param link
#'
#' @return
#' @export
#'
#' @examples
inv_link_lookup <- function(link, family = NULL) {
  if(!is.null(family)){
    switch (family,
            "logitnormal" = return(bayesfam::inv_logit),
            "cauchitnormal" = return(bayesfam::inv_cauchit),
            "cloglognormal" = return(bayesfam::inv_cloglog),
            "lognormal" = return(exp),
            "softplusnormal" = return(bayesfam::inv_softplus)
    )
  }
  switch(link,
    "logit" = bayesfam::inv_logit,
    "cauchit" = bayesfam::inv_cauchit,
    "cloglog" = bayesfam::inv_cloglog,
    "identity" = identity,
    "log" = exp,
    "softplus" = bayesfam::inv_softplus
  )
}

#' Title
#'
#' @param family
#'
#' @return
#' @export
#'
#' @examples
second_family_parameter_lookup <- function(family) {
  brms_family_lookup(family)$dpars[2]
}

#' Title
#'
#' @param family
#'
#' @return
#' @export
#'
#' @examples
prior_lookup <- function(family) {
  switch(family,
    "frechet" = c(
      brms::set_prior("", class = "Intercept"),
      brms::set_prior("", class = "nu", lb = 1.00001)
    ),
    c(
      brms::set_prior("", class = "Intercept"),
      brms::set_prior("", class = second_family_parameter_lookup(family))
    )
  )
}

#' Generate lookup keys for fit configurations to retrieve prefit objects
#' matching said config.
#'
#' @param fit_conf
#'
#' @return A hash generated from the fit configuration
#' @export
#'
#' @examples
#' fit_conf_key(
#'   list(
#'     fit_family = "gaussian",
#'     fit_link = "identity",
#'     prior = list(c(brms::set_prior("", class = "Intercept")))
#'   )
#' )
fit_conf_key <- function(fit_conf) {
  return(
    rlang::hash(
      list(
        fit_conf$fit_family,
        fit_conf$fit_link,
        fit_conf$prior
      )
    )
  )
}

#' Lookup for limits of family auxiliary parameters.
#'
#' @param family The identifier string of a family.
#'
#' @return List containing lower and upper bounds for the auxiliary parameter.
#' @export
#'
#' @examples
aux_limits_lookup <- function(family) {
  switch(family,
    "beta" = list(lb = 0, ub = Inf),
    "kumaraswamy" = list(lb = 0, ub = Inf),
    "logitnormal" = list(lb = 0, ub = Inf),
    "cauchitnormal" = list(lb = 0, ub = Inf),
    "cloglognormal" = list(lb = 0, ub = Inf),
    "simplex" = list(lb = -Inf, ub = Inf),
    "gaussian" = list(lb = 0, ub = Inf),
    "gamma" = list(lb = 0, ub = Inf),
    "weibull" = list(lb = 0, ub = Inf),
    "lognormal" = list(lb = 0, ub = Inf),
    "softplusnormal" = list(lb = 0, ub = Inf),
    "lomax" = list(lb = 1, ub = Inf),
    "frechet" = list(lb = 1, ub = Inf),
    "inverse.gaussian" = list(lb = 0, ub = Inf),
    "betaprime" = list(lb = 0, ub = Inf),
    "gompertz" = list(lb = 0, ub = Inf)
  )
}
