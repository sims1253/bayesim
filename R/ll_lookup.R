#' Title
#'
#' @param family
#'
#' @return
#' @export
#'
#' @examples
prior_lookup <- function(family) {
  switch(family,
    "frechet" = c(
      brms::set_prior("", class = "Intercept"),
      brms::set_prior("", class = "nu", lb = 1.00001)
    ),
    c(
      brms::set_prior("", class = "Intercept"),
      lapply(
        bayesfam::aux_family_parameters_lookup(family),
        function(x) brms::set_prior("", class = x)
      )
    )
  )
}

#' Generate lookup keys for fit configurations to retrieve prefit objects
#' matching said config.
#'
#' @param fit_conf
#'
#' @return A hash generated from the fit configuration
#' @export
#'
#' @examples
#' fit_conf_key(
#'   list(
#'     fit_family = "gaussian",
#'     fit_link = "identity",
#'     prior = list(c(brms::set_prior("", class = "Intercept")))
#'   )
#' )
fit_conf_key <- function(fit_conf) {
  return(
    rlang::hash(
      list(
        fit_conf$fit_family,
        fit_conf$fit_link,
        fit_conf$formula,
        fit_conf$prior
      )
    )
  )
}
