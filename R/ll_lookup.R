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
    "beta" = return(brms::Beta(link = link)),
    "kumaraswamy" = return(kumaraswamy(link = link)),
    "logitnormal" = return(logitnormal(link = link)),
    "cauchitnormal" = return(cauchitnormal(link = link)),
    "cloglognormal" = return(cloglognormal(link = link)),
    "simplex" = return(simplex(link = link)),
    "gaussian" = return(gaussian(link = link))
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
    "beta" = return(rbeta_custom),
    "kumaraswamy" = return(rkumaraswamy),
    "logitnormal" = return(rlogitnormal),
    "cauchitnormal" = return(rcauchitnormal),
    "cloglognormal" = return(rcloglognormal),
    "simplex" = return(rsimplex),
    "gaussian" = return(rnorm)
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
    "logit" = return(inv_logit),
    "cauchit" = return(inv_cauchit),
    "cloglog" = return(inv_cloglog),
    "identity" = return(identity)
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
    "beta" = return("phi"),
    "kumaraswamy" = return("p"),
    "logitnormal" = return("sigma"),
    "cauchitnormal" = return("sigma"),
    "cloglognormal" = return("sigma"),
    "simplex" = return("sigma"),
    "gaussian" = return("sigma")
  )
}
