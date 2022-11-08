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
inv_link_lookup <- function(link) {
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
