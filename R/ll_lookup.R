#' Lookup default priors for a model family
#'
#' @param family A string specifying the model family (e.g., "gaussian", "student_t")
#'
#' @return A list of brms prior objects for the specified family
#' @export
#' @keywords internal
#'
prior_lookup <- function(family) {
  aux_params <- tryCatch(
    bayesfam::aux_family_parameters_lookup(family),
    error = function(e) {
      # Fallback to common auxiliary parameters
      switch(
        family,
        "gaussian" = "sigma",
        "student_t" = c("sigma", "nu"),
        "gamma" = "shape",
        "beta" = c("phi", "beta"),
        character(0)
      )
    }
  )

  switch(
    family,
    "frechet" = c(
      brms::set_prior("", class = "Intercept"),
      brms::set_prior("", class = "nu", lb = 1.00001)
    ),
    c(
      brms::set_prior("", class = "Intercept"),
      lapply(
        aux_params,
        function(x) brms::set_prior("", class = x)
      )
    )
  )
}

#' Generate lookup keys for fit configurations to retrieve prefit objects
#' matching said config.
#'
#' @param fit_conf A list containing fit configuration with elements:
#'   fit_family, fit_link, formula, and prior
#'
#' @return A hash generated from the fit configuration
#' @export
#' @keywords internal
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
  rlang::hash(
    list(
      fit_conf$fit_family,
      fit_conf$fit_link,
      fit_conf$formula,
      fit_conf$prior
    )
  )
}
