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
  p_mean <- mean(as.vector(as_draws_array(fit, variable = "b_x")))
  return(p_mean - x_y_coef)
}
