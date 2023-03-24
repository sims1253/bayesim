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
    "gompertz" = bayesfam::gompertz(link = link),
    "student" = brms::brmsfamily("student", link = link),
    "skew_normal" = brms::brmsfamily("skew_normal", link = link),
    "generalized_normal" = bayesfam::generalized_normal(link = link),
    "asym_laplace" = brms::brmsfamily("asym_laplace", link = link),
    "exgaussian" = brms::brmsfamily("exgaussian", link = link),
    "gumbel" = bayesfam::gumbel_mean(link = link),
    "symlognormal" = bayesfam::symlognormal(link = link),
    "logistic" = bayesfam::logistic(link = link)
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
    "gompertz" = bayesfam::rgompertz,
    "logistic" = bayesfam::rlogistic,
    "student" = bayesfam::rstudent_custom,
    "skew_normal" = brms::rskew_normal,
    "generalized_normal" = bayesfam::rgeneralized_normal,
    "asym_laplace" = brms::rasym_laplace,
    "exgaussian" = bayesfam::rexgauss_custom,
    "gumbel" = bayesfam::rgumbel_mean,
    "symlognormal" = bayesfam::rsymlognormal
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
link_lookup <- function(link, family = NULL, inv = FALSE) {
  if (inv) {
    if (!is.null(family)) {
      switch(family,
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
  } else {
    switch(family,
      "logitnormal" = return(bayesfam::logit),
      "cauchitnormal" = return(bayesfam::cauchit),
      "cloglognormal" = return(bayesfam::cloglog),
      "lognormal" = return(log),
      "softplusnormal" = return(bayesfam::softplus)
    )
    switch(link,
      "logit" = bayesfam::logit,
      "cauchit" = bayesfam::cauchit,
      "cloglog" = bayesfam::cloglog,
      "identity" = identity,
      "log" = log,
      "softplus" = bayesfam::softplus
    )
  }
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
    # maybe implement this more generalized, to facilitate other custom
    # families with 3 or more parameters
    "generalized_normal" = c(
      brms::set_prior("", class = "Intercept"),
      brms::set_prior("", class = "sigma", lb = 0),
      brms::set_prior("", class = "beta", lb = 0)
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
        fit_conf$formula,
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
    "gompertz" = list(lb = 0, ub = Inf),
    "student" = list(lb = c(0, 0), ub = c(Inf, Inf)),
    "skew_normal" = list(lb = c(0, -Inf), ub = c(Inf, Inf)),
    "generalized_normal" = list(lb = c(0, 0), ub = c(Inf, Inf)),
    "asym_laplace" = list(lb = c(0, 0), ub = c(Inf, 1)),
    "exgaussian" = list(lb = c(0, 0), ub = c(Inf, Inf)),
    "gumbel" = list(lb = 0, ub = Inf),
    "symlognormal" = list(lb = 0, ub = Inf),
    "logistic" = list(lb = 0, ub = Inf)
  )
}
