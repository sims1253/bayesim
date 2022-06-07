#' Title
#'
#' @param family
#' @param link
#'
#' @return
#' @export
#'
#' @examples
brms_family_lookup <- function(family, link) {
  switch(family,
    "beta" = brms::brmsfamily("beta", link = link),
    "kumaraswamy" = kumaraswamy(link = link),
    "logitnormal" = logitnormal(link = link),
    "cauchitnormal" = cauchitnormal(link = link),
    "cloglognormal" = cloglognormal(link = link),
    "simplex" = simplex(link = link),
    "gaussian" = brms::brmsfamily("gaussian", link = link),
    "gamma" = brms::brmsfamily("gamma", link = link),
    "weibull" = brms::brmsfamily("weibull", link = link),
    "lognormal" = brms::brmsfamily("lognormal", link = link),
    "softplusnormal" = softplusnormal(link = link),
    "lomax" = lomax(link = link),
    "frechet" = brms::brmsfamily("frechet", link = link),
    "inverse.gaussian" = brms::brmsfamily("inverse.gaussian", link = link),
    "betaprime" = betaprime(link = link),
    "gompertz" = gompertz(link = link)
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
rng_lookup <- function(family, link = NULL) {
  switch(family,
    "beta" = rbeta_custom,
    "kumaraswamy" = rkumaraswamy,
    "logitnormal" = rlogitnormal,
    "cauchitnormal" = rcauchitnormal,
    "cloglognormal" = rcloglognormal,
    "simplex" = rsimplex,
    "gaussian" = rnorm,
    "gamma" = rgamma_custom,
    "weibull" = rweibull_custom,
    "lognormal" = rlognormal,
    "softplusnormal" = rsoftplusnormal,
    "lomax" = rlomax,
    "frechet" = rfrechet_custom,
    "inverse.gaussian" = brms::rinv_gaussian,
    "betaprime" = rbetaprime,
    "gompertz" = rgompertz
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
inv_link_lookup <- function(link) {
  switch(link,
    "logit" = inv_logit,
    "cauchit" = inv_cauchit,
    "cloglog" = inv_cloglog,
    "identity" = identity,
    "log" = exp,
    "softplus" = inv_softplus
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
  switch(family,
    "beta" = "phi",
    "kumaraswamy" = "p",
    "logitnormal" = "sigma",
    "cauchitnormal" = "sigma",
    "cloglognormal" = "sigma",
    "simplex" = "sigma",
    "gaussian" = "sigma",
    "gamma" = "shape",
    "weibull" = "shape",
    "lognormal" = "sigma",
    "softplusnormal" = "sigma",
    "lomax" = "alpha",
    "frechet" = "nu",
    "inverse.gaussian" = "shape",
    "betaprime" = "phi",
    "gompertz" = "eta"
  )
}
