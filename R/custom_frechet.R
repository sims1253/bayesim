#' Title
#'
#' @param x
#' @param mu Mean
#' @param k Shape
#'
#' @return
#' @export
#'
#' @examples
dfrechet_custom <- function(x, mu, a) {
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
  return(brms::dfrechet(x = x, loc = 0, shape = k, scale = mu / gamma(1 + 1 / k)))
}


#' Title
#'
#' @param n
#' @param mu Mean
#' @param k Shape
#'
#' @return
#' @export
#'
#' @examples
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
