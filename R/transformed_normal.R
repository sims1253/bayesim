#' Title
#'
#' @param link
#' @param link_sigma
#'
#' @return
#' @export
#'
#' @examples
transformed_normal <- function(link = "logit", link_sigma = "log") {
  if (link == "logit") {
    return(logitnormal(link_sigma = link_sigma))
  }
  if (link == "cauchit") {
    return(cauchitnormal(link_sigma = link_sigma))
  }
  if (link == "cloglog") {
    return(cloglognormal(link_sigma = link_sigma))
  }
  stop(paste(link, "is not a valid link for the transformed_normal family."))
}
