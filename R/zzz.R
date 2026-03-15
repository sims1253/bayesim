#' @keywords internal
#' @importFrom utils globalVariables
utils::globalVariables(c("sim_id", "par_seed", "%>%"))

.onLoad <- function(libname, pkgname) {
  # Register built-in metrics on package load
  register_built_in_metrics()
}
