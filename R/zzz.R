.onLoad <- function(libname, pkgname) {
  # Register built-in metrics on package load
  register_built_in_metrics()
}
