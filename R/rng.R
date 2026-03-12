#' Set up global RNG for simulation
#'
#' Sets RNG kind to L'Ecuyer-CMRG and initializes with global seed.
#' Must be called once at simulation start.
#'
#' @param seed Integer seed for the simulation
#'
#' @return The initial `.Random.seed` state (invisibly)
#'
#' @export
#'
#' @examples
#' setup_global_rng(42)
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
#' @export
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

#' Advance RNG stream
#'
#' Advances the RNG state by n steps without returning random values.
#' Useful for skipping ahead in a stream.
#'
#' @param rng_stream Integer vector RNG state
#' @param n Number of steps to advance (default: 1)
#'
#' @return Advanced RNG state (integer vector)
#'
#' @details
#' This is a **pure function** with no side effects:
#' - It does NOT modify `.Random.seed` in the global environment
#' - It creates a local copy of the RNG state, advances it, and returns the new state
#' - The caller is responsible for setting the returned state if needed
#'
#' @export
#'
#' @examples
#' \dontrun{
#' streams <- create_task_rng_streams(42, 10)
#' advanced <- advance_rng_stream(streams[[1]], n = 5)
#' }
advance_rng_stream <- function(rng_stream, n = 1L) {
  .Random.seed <- rng_stream
  for (i in seq_len(n)) {
    runif(1)
  }
  get(".Random.seed", inherits = FALSE)
}
