# brms_model() and model_grid() — ergonomics for building brms fit_grid rows.
# D3: hand-building list-columns (fit_grid$formula <- list(...)) is the ugliest
# UX in the package; these helpers produce a tidy tibble of model specs.

#' Construct a single brms model specification
#'
#' Builds a one-row list of brms model components (formula, family, prior,
#' and stanvars suitable for assembling a fit_grid via [model_grid()].
#' Inputs are validated at construction so errors surface before compilation.
#'
#' @param formula A formula or brmsformula.
#' @param family A brms family (e.g. `gaussian()`, `student()`), or NULL for the
#'   brms default.
#' @param prior A brms prior object, or NULL.
#' @param stanvars A brms `stanvars` object, or NULL.
#'
#' @return A named list with elements `formula`, `family`, `prior`, `stanvars`.
#' @export
#' @seealso [model_grid()]
#' @examples
#' \dontrun{
#' brms_model(y ~ x, family = gaussian())
#' }
brms_model <- function(
  formula,
  family = NULL,
  prior = NULL,
  stanvars = NULL
) {
  if (
    is.null(formula) ||
      (!inherits(formula, "formula") &&
        !inherits(formula, "brmsformula"))
  ) {
    stop(bayesim_config_error(
      "formula must be a formula or brmsformula"
    ))
  }
  list(
    formula = formula,
    family = family,
    prior = prior,
    stanvars = stanvars
  )
}

#' Assemble a tibble of model specifications for a fit_grid
#'
#' Combines named [brms_model()] specs into a single data frame with one row per
#' model and `formula`/`family`/`prior`/`stanvars` list-columns,
#' ready to pass as `fit_grid` to [simulation_config()]. A `model` column holds
#' the spec names and lands in the summary as `fit_model`.
#'
#' @param ... Named `brms_model()` specs (or any named lists with the same
#'   component names).
#'
#' @return A tibble with columns `model`, `formula`, `family`, `prior`,
#'   `stanvars`.
#' @export
#' @seealso [brms_model()], [simulation_config()]
#' @examples
#' \dontrun{
#' grid <- model_grid(
#'   gaussian = brms_model(y ~ x, gaussian()),
#'   student  = brms_model(y ~ x, student())
#' )
#' }
model_grid <- function(...) {
  specs <- list(...)
  if (length(specs) == 0L) {
    stop(bayesim_config_error("model_grid() requires at least one named spec"))
  }
  if (is.null(names(specs)) || any(names(specs) == "")) {
    stop(bayesim_config_error("every model_grid() spec must be named"))
  }
  components <- c("formula", "family", "prior", "stanvars")
  for (nm in names(specs)) {
    spec <- specs[[nm]]
    if (!is.list(spec) || !all(components %in% names(spec))) {
      stop(bayesim_config_error(paste(
        "spec",
        nm,
        "must be a brms_model() (or a list with components:",
        paste(components, collapse = ", "),
        ")"
      )))
    }
  }

  out <- tibble::tibble(
    model = names(specs),
    formula = lapply(specs, `[[`, "formula"),
    family = lapply(specs, `[[`, "family"),
    prior = lapply(specs, `[[`, "prior"),
    stanvars = lapply(specs, `[[`, "stanvars")
  )
  out
}
