#' Title
#'
#' @param family
#'
#' @return
#' @export
#'
#' @examples
get_prefit <- function(family) {
  prefit <- brms::brm(
    y ~ 1 + x,
    data = list(y = c(0.5), x = c(1)),
    family = family,
    stanvars = family$stanvars,
    chains = 0,
    refresh = 0,
    silent = 2,
    control = list(adapt_delta = 0.9),
    prior = c(
      brms::prior("", class = "Intercept")
    )
  )
  return(prefit)
}
