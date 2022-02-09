#' Title
#'
#' @param fit
#' @param x_y_coef
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
p_bias <- function(fit, x_y_coef, ...) {
  return(p_mean(fit) - x_y_coef)
}

#' Title
#'
#' @param fit
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
p_mean <- function(fit) {
  return(mean(as.vector(as_draws_array(fit, variable = "b_x"))))
}

#' Title
#'
#' @param fit
#'
#' @return
#' @export
#'
#' @examples
p_rstar <- function(fit) {
  posterior::rstar(as_draws_array(fit))
}
