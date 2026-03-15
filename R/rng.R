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
#' @keywords internal
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
#' Advances the RNG state by n steps and returns the new state.
#'
#' @param rng_stream Integer vector RNG state
#' @param n Number of steps to advance (default: 1)
#'
#' @return Advanced RNG state (integer vector)
#'
#' @details
#' This function temporarily sets `.Random.seed` in `.GlobalEnv`
#' to advance the stream, then returns the new state. It does NOT permanently
#' modify the global `.Random.seed` (unless it did not exist before, in which
#' case it is removed after use). The caller is responsible for applying
#' the returned state if needed (e.g., via `set_task_rng()`).
#'
#' @keywords internal
#' @export
#'
#' @examples
#' \dontrun{
#' streams <- create_task_rng_streams(42, 10)
#' advanced <- advance_rng_stream(streams[[1]], n = 5)
#' }
advance_rng_stream <- function(rng_stream, n = 1L) {
  old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  on.exit(
    if (seed_existed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv, inherits = FALSE)
    } else {
      rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    },
    add = TRUE
  )
  assign(".Random.seed", rng_stream, envir = .GlobalEnv, inherits = FALSE)
  for (i in seq_len(n)) {
    runif(1)
  }
  get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
}
