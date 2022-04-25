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
rbeta_custom <- function(n, mu, phi) {
  if (isTRUE(any(mu <= 0 || mu >= 1))) {
    stop("The median must be in (0,1).")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("P must be above 0.")
  }
  return(rbeta(n, mu * phi, (1 - mu) * phi))
}

#' Title
#'
#' @param x
#' @param mu
#' @param phi
#'
#' @return
#' @export
#'
#' @examples
dbeta_custom <- function(x, mu, phi) {
  dbeta(x, mu * phi, (1 - mu) * phi)
}
