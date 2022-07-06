#' Frechet density function in mean parametrisation.
#'
#' @param x
#' @param mu Mean
#' @param nu Shape
#'
#' @details Define scale parameter sigma as
#' \deqn{\sigma(\mu, k) := \mu / \Gamma(1 + 1 / k)}
#' @details The Frechet distribution has density
#' \deqn{f(y) = (k /\sigma) * (y / \sigma)^{-(1 + k)} * exp(-(x / \sigma)^{-k}) }
#'
#' @return dfrechet(x | mu, k)
#' @export
#'
#' @examples
dfrechet_custom <- function(x, mu, nu) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("frechet is only defined for x > 0")
  }
  if (isTRUE(nu <= 1)) {
    stop("frechet is only defined for nu > 1")
  }
  if (isTRUE(mu <= 0)) {
    stop("frechet is only defined for mu > 0")
  }
  return(brms::dfrechet(x = x, loc = 0, scale = mu / gamma(1 - 1 / nu), shape = nu))
}


#' Frechet RNG function in Mean parametrisation
#'
#' @param n
#' @param mu Mean
#' @param nu Shape
#'
#' @return n samples as Frechet-distribution
#' @export
#'
#' @examples
rfrechet_custom <- function(n, mu, nu) {
  # check the arguments
  if (isTRUE(nu <= 1)) {
    stop("frechet is only defined for nu > 1")
  }
  if (isTRUE(mu <= 0)) {
    stop("frechet is only defined for mu > 0")
  }
  return(brms::rfrechet(n = n, scale = mu / gamma(1 - 1 / nu), shape = nu))
}
