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
    return(beta_custom(link = link))
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
}
