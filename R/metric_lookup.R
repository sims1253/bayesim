#' Title
#'
#' @param identifier
#' @param fit
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
metric_lookup <- function(identifier, fit, ...) {
  if (grepl("^hyp_", identifier)) {
    hyp <- substring(identifier, 5)
    return(hypothesis_cis(hyp, fit, ...))
  }

  if (identifier == "rmse_s") {
    return(rmse_s(fit, ...))
  }

  if (identifier == "bias") {
    return(p_bias(fit, ...))
  }
}
