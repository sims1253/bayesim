#' RNG for the beta distribution using brms' mean and phi parametrization.
#'
#' @param n
#' @param mu
#' @param phi
#'
#' @return
#' @export
#'
#' @examples
rBeta <- function(n, mu, phi) {
  rbeta(n, mu * phi, (1 - mu) * phi)
}
