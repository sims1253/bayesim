#' Custom Gamma density function in mean parametrisation.
#'
#' @param x Value space, x > 0.
#' @param mu Mean parameter of the density, mu > 0.
#' @param k Shape parameter, k > 0.
#'
#' @details Define rate constante rho as:
#' \deqn{\rho(\alpha, \mu) = \frac{\alpha}{\mu}}
#' @details The Frechet distribution density is defined as
#' \deqn{f(y) = \frac{y^{\alpha - 1} exp(-\frac{y}{\rho})} {\rho^\alpha \Gamma(\alpha)} }
#'
#' @return f(x | mu, k)
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 10, length.out = 1000)
#' plot(x, dgamma_custom(x, mu = 2, a = 2), type = "l")
dgamma_custom <- function(x, mu, a, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("gamma is only defined for x > 0")
  }
  if (isTRUE(a <= 0)) {
    stop("gamma is only defined for a > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("gamma is only defined for mu > 0")
  }
  return(dgamma(x = x, shape = a, rate = a / mu, log = log))
}


#' Custom Gamma RNG function in mean parametrisation.
#'
#' @param n Number of samples, scalar natural number.
#' @param mu Mean parameter, mu > 0.
#' @param k Shape parameter, k > 0.
#'
#' @return n Gamma distributed samples.
#' @export
#'
#' @examples hist(log(rgamma_custom(10000, mu = 2, a = 1)))
rgamma_custom <- function(n, mu, a) {
  # check the arguments
  if (isTRUE(a <= 0)) {
    stop("gamma is only defined for a > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("gamma is only defined for mu > 0")
  }
  return(rgamma(n = n, shape = a, rate = a / mu))
}
