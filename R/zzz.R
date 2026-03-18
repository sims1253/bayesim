.onLoad <- function(libname, pkgname) {
  # Register built-in metrics on package load
  register_built_in_metrics()
}

.onUnload <- function(libpath) {
  # Clear the metric registry to avoid polluting the R session
  clear_registry()
}
