#' Title
#'
#' @param x
#' @param mu
#' @param alpha
#' @param log
#'
#' @return
#' @export
#'
#' @examples
dlomax <- function(x, mu, alpha, log = FALSE) {
  # check arguments
  if(isTRUE(any(x < 0))) {
    stop("The Lomax PDF is defined only on the positive scale")
  }
  if(isTRUE(alpha <= 1)) {
    stop("The Lomax PDF is defined only for alpha > 1")
  }
  if(isTRUE(mu <= 0)) {
    stop("The Lomax PDF is defined only for mu > 0")
  }

  lpdf <- log(alpha) +
          alpha * (log(mu) + log(alpha - 1)) -
          (alpha + 1) * (log(x + (mu * (alpha - 1))))

  if(log) {
    return(lpdf)
  }
  else {
    return(exp(lpdf))
  }
}
