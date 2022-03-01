#' Title
#'
#' @param identifier
#' @param link
#'
#' @return
#' @export
#'
#' @examples
likelihood_lookup <- function(identifier, link) {
  if (identifier == "beta") {
    return(brms::Beta(link = link))
  }

  if (identifier == "kumaraswamy") {
    return(kumaraswamy(link = link))
  }

  if (identifier == "logitnormal") {
    return(logitnormal(link = link))
  }

  if (identifier == "simplex") {
    return(simplex(link = link))
  }

  if (identifier == "gaussian") {
    return(gaussian(link = link))
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
rng_lookup <- function(identifier) {
  if (identifier == "beta") {
    return(rbeta_custom)
  }

  if (identifier == "kumaraswamy") {
    return(rkumaraswamy)
  }

  if (identifier == "logitnormal") {
    return(rlogitnormal)
  }

  if (identifier == "simplex") {
    return(rsimplex)
  }

  if (identifier == "gaussian") {
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
    return(logistic)
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
