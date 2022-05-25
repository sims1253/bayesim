#' Logit link function
#'
#' @param x value of x to be transformed, 0 < x < 1
#'
#' @return logit value of x
#' @export
#'
#' @examples x <- seq(from = 0.1 , to = 0.9 , length.out = 100)
#' y <- logit(x)
logit <- function(x) {
  return(qlogis(x))
}

#' Inverse-Logit link function, equal to Logistic link
#'
#' @param x value of x to be transformed, any real scalar or vector allowed
#'
#' @return inverse logit value of x
#' @export
#'
#' @examples x <- seq(from = -100 , to = 100 , length.out = 1000)
#' y <- inv_logit(x)
inv_logit <- function(x) {
  return(plogis(x))
}

#' Logistic link function, equivalent to Inverse-Logit link
#'
#' @param x value of x to be transformed, any real scalar or vector allowed
#'
#' @return logistic value of x
#' @export
#'
#' @examples x <- seq(from = -100 , to = 100 , length.out = 1000)
#' y <- logistic(x)
logistic <- function(x) {
  return(inv_logit(x))
}


#' Complementary-Log-Log link function
#'
#' @param x value of x to be transformed, 0 < x < 1
#'
#' @return cloglog value of x
#' @export
#'
#' @examples x <- seq(from = 0.1 , to = 0.9 , length.out = 100)
#' y <- cloglog(x)
cloglog <- function(x) {
  log(-log1p(-x))
}


#' Inverse CLogLog link function
#'
#' @param x value of x to be transformed, any real scalar or vector allowed
#'
#' @return inverse-cloglog value of x
#' @export
#'
#' @examples x <- seq(from = -3, to = 1 , length.out = 100)
#' y <- inv_cloglog(x)
inv_cloglog <- function(x) {
  return(1 - exp(-exp(x)))
}


#' Cauchit link function
#'
#' @param x value of x to be transformed, not defined for x of whole integer
#'
#' @return cauchit value of x
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 0.9 , length.out = 100)
#' y <- cauchit(x)
cauchit <- function(x) {
  return(qcauchy(x))
}


#' Inverse Cauchit link function equivalent to pcauchy
#'
#' @param x value of x to be transformed, any real scalar or vector allowed
#'
#' @return inverse cauchit of x
#' @export
#'
#' @examples x <- seq(from = -10 , to = 10 , length.out = 100)
#' y <- inv_cauchit(x)
inv_cauchit <- function(x) {
  return(pcauchy(x))
}
