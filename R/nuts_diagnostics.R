#' Title
#'
#' @param fit
#'
#' @return
#' @export
#'
#' @examples
divergents <- function(fit) {
  nuts_diag <- brms::nuts_params(fit)
  return(sum(nuts_diag$Value[which(nuts_diag$Parameter == "divergent__")]))
}
