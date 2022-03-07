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
  if (family == "beta") {
    return(brms::Beta(link = link))
  }

  if (family == "kumaraswamy") {
    return(kumaraswamy(link = link))
  }

  if (family == "transformed_normal") {
    return(transformed_normal(link = link))
  }

  if (family == "simplex") {
    return(simplex(link = link))
  }

  if (family == "gaussian") {
    return(gaussian(link = link))
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
rng_lookup <- function(family, link = NULL) {
  if (family == "beta") {
    return(rbeta_custom)
  }

  if (family == "kumaraswamy") {
    return(rkumaraswamy)
  }

  if (family == "transformed_normal") {
    switch(link,
      "logit" = return(rlogitnormal),
      "cauchit" = return(rcauchitnormal),
      "cloglog" = return(rcloglognormal)
    )
  }

  if (family == "simplex") {
    return(rsimplex)
  }

  if (family == "gaussian") {
    return(rnorm)
  }
}


#' Title
#'
#' @param identifier
#'
#' @return
#' @export
#'
#' @examples
inv_link_lookup <- function(link) {
  if (link == "logit") {
    return(inv_logit())
  }

  if (link == "cloglog") {
    return(inv_cloglog)
  }

  if (link == "cauchit") {
    return(inv_cauchit)
  }

  if (link == "identity") {
    return(function(x) {
      x
    })
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
  if (family == "beta") {
    return("phi")
  }

  if (family == "kumaraswamy") {
    return("p")
  }

  if (family == "transformed_normal") {
    return("sigma")
  }

  if (family == "simplex") {
    return("sigma")
  }

  if (family == "gaussian") {
    return("sigma")
  }
}
