# Test Fixtures and Helpers for bayesim Tests
#
# This file provides reusable test fixtures for testing bayesim components.
# These helpers create consistent, valid test data structures.

# Null Coalescing Operator ------------------------------------------------

#' Null coalescing operator
#'
#' Returns the left-hand side if not NULL, otherwise the right-hand side.
#' Useful for providing default values in test fixtures.
#'
#' @param x Value to check for NULL
#' @param y Default value to return if x is NULL
#' @return x if not NULL, otherwise y
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

# Mock Data Generator -----------------------------------------------------

#' Mock data generator for testing
#'
#' Creates synthetic training and test data for simulation testing.
#' Generates simple linear regression data with known parameters.
#'
#' @param data_spec A list specifying data parameters (e.g., n for sample size)
#' @param seed Integer seed for reproducible random data generation
#' @param task_ctx Task context containing task_id, data_idx, fit_idx, rep_idx
#'
#' @return A list containing:
#'   \describe{
#'     \item{train}{Data frame with y and x columns for training}
#'     \item{test}{NULL (no test data in mock)}
#'     \item{response}{Name of the response variable}
#'     \item{true_params}{Named vector of true parameter values}
#'     \item{vars_of_interest}{Names of parameters to estimate}
#'     \item{meta}{List of metadata}
#'   }
#'
#' @examples
#' mock_data_generator(list(n = 50), list(task_id = "test"))
#'
#' @keywords internal
mock_data_generator <- function(data_spec, task_ctx) {
  # Consume the ambient RNG state (the worker restores the per-task L'Ecuyer
  # stream before each call); do not re-seed internally.
  n <- data_spec$n %||% 100
  list(
    train = data.frame(
      y = stats::rnorm(n),
      x = stats::rnorm(n)
    ),
    test = NULL,
    response = "y",
    true_params = c(beta = 1.0, sigma = 1.0),
    vars_of_interest = c("beta", "sigma"),
    meta = list()
  )
}

# Valid Test Fixtures ------------------------------------------------------

#' Create a valid data bundle for testing
#'
#' Generates a minimal valid data bundle that satisfies the data contract.
#' Useful for testing functions that accept data bundles as input.
#'
#' @return A list containing a valid data bundle with:
#'   \describe{
#'     \item{train}{Data frame with y and x columns (10 rows)}
#'     \item{test}{NULL (no test data)}
#'     \item{response}{Name of the response variable ("y")}
#'     \item{true_params}{Named vector of true parameter values}
#'     \item{vars_of_interest}{Names of parameters to estimate}
#'     \item{meta}{List of metadata}
#'   }
#'
#' @examples
#' bundle <- valid_data_bundle()
#' stopifnot(!is.null(bundle$train))
#' stopifnot(bundle$response == "y")
#'
#' @keywords internal
valid_data_bundle <- function() {
  list(
    train = data.frame(y = 1:10, x = 1:10),
    test = NULL,
    response = "y",
    true_params = c(beta = 1.0),
    vars_of_interest = "beta",
    meta = list()
  )
}

#' Create a valid fit result for testing
#'
#' Generates a valid fit result object with simulated posterior draws.
#' Useful for testing functions that process or validate fit results.
#'
#' @return A fit result object (created via new_fit_result) containing:
#'   \describe{
#'     \item{success}{TRUE indicating successful fit}
#'     \item{draws}{50x2 matrix of posterior draws (alpha, beta)}
#'     \item{diagnostics}{List with convergence metrics (rhat = 1.01)}
#'     \item{timing}{List with timing info (total, warmup, sample)}
#'   }
#'
#' @examples
#' fit <- valid_fit_result()
#' stopifnot(fit$success)
#' stopifnot(ncol(fit$draws) == 2)
#'
#' @keywords internal
valid_fit_result <- function() {
  draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
  colnames(draws) <- c("alpha", "beta")
  new_fit_result(
    success = TRUE,
    draws = draws,
    diagnostics = list(rhat = 1.01),
    timing = list(total = 1.0, warmup = 0.5, sample = 0.5)
  )
}

#' Create a valid task context for testing
#'
#' Generates a valid task context with standard indices.
#' Useful for testing functions that require task context information.
#'
#' @return A list containing:
#'   \describe{
#'     \item{task_id}{Unique task identifier string}
#'     \item{data_idx}{Data configuration index (integer)}
#'     \item{fit_idx}{Fit configuration index (integer)}
#'     \item{rep_idx}{Replication index (integer)}
#'   }
#'
#' @examples
#' ctx <- valid_task_ctx()
#' stopifnot(ctx$task_id == "d001_f001_r00001")
#' stopifnot(ctx$data_idx == 1)
#'
#' @keywords internal
valid_task_ctx <- function() {
  list(
    task_id = "d001_f001_r00001",
    data_idx = 1,
    fit_idx = 1,
    rep_idx = 1
  )
}
