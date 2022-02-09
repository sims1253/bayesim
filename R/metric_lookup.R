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

  if (identifier == "divergents") {
    return(divergents(fit))
  }

  if (identifier == "pmean") {
    return(p_mean(fit))
  }

  if (identifier == "rstar") {
    return(p_rstar(fit))
  }
}
