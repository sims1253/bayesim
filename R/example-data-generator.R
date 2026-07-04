# A minimal built-in data generator.
#
# This is the simplest package-internal generator: a fixed-truth linear
# regression. It is useful as a default/reference pattern, for documentation
# examples, and for the manifest-rehydration tests (because it lives in the
# package namespace and can be rehydrated from a checkpoint manifest).

#' Built-in example data generator
#'
#' Generates a simple linear-regression dataset (`y = beta * x + noise`).
#'
#' Generators consume the ambient RNG state (the worker restores the per-task
#' L'Ecuyer stream via `set_task_rng()`); the `seed` argument is retained in
#' the signature for backends requiring an explicit integer but is not used to
#' re-seed.
#'
#' @param data_spec List with `n`, `beta`, `sigma`.
#' @param seed Scalar task seed (unused; ambient RNG is consumed).
#' @param task_ctx Task context.
#' @return A `data_bundle` list.
#' @keywords internal
bayesim_example_data_generator <- function(data_spec, seed, task_ctx) {
  n <- as.integer(data_spec$n %||% 10L)
  beta <- as.numeric(data_spec$beta %||% 1)
  sigma <- as.numeric(data_spec$sigma %||% 1)

  x <- stats::rnorm(n)
  y <- beta * x + stats::rnorm(n, sd = sigma)

  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(beta = beta, sigma = sigma),
    vars_of_interest = c("beta", "sigma"),
    references = c(beta = 0, sigma = 1),
    meta = list(task_id = task_ctx$task_id %||% NA_character_)
  )
}
