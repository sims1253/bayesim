.onLoad <- function(libname, pkgname) {
  # bayesim 2.0 uses constructors directly; no metric registry.
  # F3: register as_tibble S3 method with tibble's namespace so
  # tibble::as_tibble(result) dispatches to our method.
  s3 <- list(
    list("tibble", "as_tibble", "bayesim_simulation_result", as_tibble.bayesim_simulation_result)
  )
  for (m in s3) {
    if (!requireNamespace(m[[1]], quietly = TRUE)) next
    registerS3method(m[[2]], m[[3]], m[[4]], envir = asNamespace(m[[1]]))
  }
  invisible()
}
