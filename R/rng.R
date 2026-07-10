#' Set up global RNG for simulation
#'
#' Sets RNG kind to L'Ecuyer-CMRG and initializes with global seed.
#' Must be called once at simulation start.
#'
#' @param seed Integer seed for the simulation
#'
#' @return The initial `.Random.seed` state (invisibly)
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' setup_global_rng(42)
#' }
setup_global_rng <- function(seed) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  invisible(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
}

#' Set RNG state for a task
#'
#' Restores .Random.seed from a precomputed stream.
#' Call this at the start of each worker task.
#'
#' @param rng_stream Integer vector from create_task_rng_streams()
#'
#' @return NULL (invisibly). Side effect: sets `.Random.seed` in global environment.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' streams <- create_task_rng_streams(42, 10)
#' set_task_rng(streams[[1]])
#' }
set_task_rng <- function(rng_stream) {
  assign(".Random.seed", rng_stream, envir = .GlobalEnv)
  invisible(NULL)
}
