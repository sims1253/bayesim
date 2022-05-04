#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
logit <- function(x) {
  return(qlogis(x))
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
inv_logit <- function(x) {
  return(plogis(x))
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
  return(inv_logit(x))
}


#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
cloglog <- function(x) {
  log(-log1p(-x))
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
cauchit <- function(x) {
  return(qcauchy(x))
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
