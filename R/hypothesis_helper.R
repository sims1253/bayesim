#' Title
#'
#' @param hyp
#' @param fit
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
hypothesis_cis <- function(hyp, fit, ...) {
  result <- brms::hypothesis(fit, hyp)
  return(list(
    "lower_ci" = result$hypothesis$CI.Lower,
    "upper_ci" = result$hypothesis$CI.Upper,
    "post_prop" = result$hypothesis$Post.Prob
  ))
}
