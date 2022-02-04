#' Title
#'
#' @param fit
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
rmse_s <- function(fit, x_y_coef, ...) {
  posterior_draws <- as.vector(as_draws_array(fit, variable = "b_x"))
  return(sqrt(mean((posterior_draws - x_y_coef)^2)))
}
