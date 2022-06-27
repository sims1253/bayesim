#' Frechet density function in mean parametrisation.
#'
#' @param x Value space, x > 0
#' @param mu Mean parameter, mu > 0
#' @param k Shape parameter, k > 0
#'
#' @details Define scale parameter sigma as
#' \deqn{\sigma(\mu, k) := \mu / \Gamma(1 + 1 / k)}
#' @details The Frechet distribution has density
#' \deqn{f(y) = (k /\sigma) * (y / \sigma)^{-(1 + k)} * exp(-(x / \sigma)^{-k}) }
#'
#' @return dfrechet(x | mu, k)
#' @export
#'
#' @examples x <- seq(from=0.01, to=10, length.out=1000)
#' plot(x, dfrechet_custom(x, mu=2, k=1), type="l")
dfrechet_custom <- function(x, mu, k, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("frechet is only defined for x > 0")
  }
  if (isTRUE(k <= 0)) {
    stop("frechet is only defined for k > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("frechet is only defined for mu > 0")
  }
  return(brms::dfrechet(x = x, loc = 0, shape = k, scale = mu / gamma(1 + 1 / k), log))
}


#' Frechet RNG function in Mean parametrisation
#'
#' @param n Number of samples as whole integer
#' @param mu Mean parameter, mu > 0
#' @param k Shape parameter, k > 0
#'
#' @return n samples as Frechet-distribution
#' @export
#'
#' @examples hist(log(rfrechet_custom(10000, mu=2, k=1)))
rfrechet_custom <- function(n, mu, k) {
  # check the arguments
  if (isTRUE(k <= 0)) {
    stop("frechet is only defined for k > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("frechet is only defined for mu > 0")
  }
  return(brms::rfrechet(n = n, shape = k, scale = mu / gamma(1 + 1 / k)))
}
