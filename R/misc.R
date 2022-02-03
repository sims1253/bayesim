#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
logit <- function(x) {
  return(log(x) - log1p(-x))
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
logistic <- function(x) {
  return(1 / (1 + exp(-x)))
}

#' Title
#'
#' @param y
#' @param mu
#'
#' @return
#' @export
#'
#' @examples
unit_deviance <- function(y, mu) {
  return(
    ((y - mu)^2) /
      (y * (1 - y) * mu^2 * (1 - mu)^2)
  )
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
inv_cloglog <- function(x) {
  return(1 - exp(-exp(x)))
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
inv_cauchit <- function(x) {
  return(pcauchy(x))
}


#' Title
#'
#' @param n
#' @param mu
#' @param phi
#'
#' @return
#' @export
#'
#' @examples
rbeta_mu <- function(n, mu, phi) {
  rbeta(n, mu * phi, (1 - mu) * phi)
}
