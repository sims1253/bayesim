#' Logit link function
#'
#' @param x value of x to be transformed, x e (0, 1)
#'
#' @return logit value of x, x unbound
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 0.9, length.out = 100)
#' plot(x, logit(x), type="l")
logit <- function(x) {
  return(qlogis(x))
}

#' Inverse-Logit link function, equal to Logistic link
#'
#' @param x value of x to be transformed, x unbound
#'
#' @return inverse logit value of x, result is e (0, 1)
#' @export
#'
#' @examples x <- seq(from = -5, to = 5, length.out = 100)
#' plot(x, inv_logit(x), type="l")
inv_logit <- function(x) {
  return(plogis(x))
}

#' Logistic link function, equivalent to Inverse-Logit link
#'
#' @param x value of x to be transformed, unbound
#'
#' @return logistic value of x, result ise (0, 1)
#' @export
#'
#' @examples x <- seq(from = -100, to = 100, length.out = 100)
#' plot(x, logistic(x), type="l")
logistic <- function(x) {
  return(inv_logit(x))
}


#' Complementary-Log-Log link function
#'
#' @param x value of x to be transformed, x e (0, 1)
#'
#' @return cloglog value of x, result is unbound
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 0.9, length.out = 100)
#' plot(x, cloglog(x), type="l")
cloglog <- function(x) {
  log(-log1p(-x))
}


#' Inverse CLogLog link function
#'
#' @param x value of x to be transformed, x unbound
#'
#' @return inverse-cloglog value of x, result is e (0, 1)
#' @export
#'
#' @examples x <- seq(from = -3, to = 1, length.out = 100)
#' plot(x, inv_cloglog(x), type="l")
inv_cloglog <- function(x) {
  return(1 - exp(-exp(x)))
}


#' Cauchit link function
#'
#' @param x value of x to be transformed, not defined for x of whole integer
#'
#' @return cauchit value of x, result is unbound
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 0.9, length.out = 100)
#' plot(x, cauchit(x), type="l")
cauchit <- function(x) {
  return(qcauchy(x))
}


#' Inverse Cauchit link function equivalent to pcauchy
#'
#' @param x value of x to be transformed, any real scalar or vector allowed
#'
#' @return inverse cauchit of x, result is e (0, 1)
#' @export
#'
#' @examples x <- seq(from = -10, to = 10, length.out = 100)
#' plot(x, inv_cauchit(x), type="l")
inv_cauchit <- function(x) {
  return(pcauchy(x))
}

#' Gaussion Error function
#'
#' @param x value to be transformed, x unbound
#'
#' @return erf function of x, result e (0, 1)
#' @export
#'
#' @examples x <- seq(from = -2, to = 2, length.out = 100)
#' plot(x, erf(x), type="l")
erf <- function(x) {
  return(2 * pnorm(x * sqrt(2)) - 1)
}


#' Softplus link function
#'
#' @param x value to be transformed, x unbound
#'
#' @return softplus function of x, result is positive unbound
#' @export
#'
#' @examples x <- seq(from = -5, to = 5, length.out = 100)
#' plot(x, softplus(x), type="l")
softplus <- function(x) {
  return(log(exp(x) + 1))
}


#' Inverse Softplus link function
#'
#' @param x value to be transformed, x is unbound
#'
#' @return inv_softplus of x, result is positive unbound
#' @export
#'
#' @examples x <- seq(from = -5, to = 5, length.out = 100)
#' plot(x, softplus(x), type="l")
inv_softplus <- function(x) {
  return(log(exp(x) - 1))
}
