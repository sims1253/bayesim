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

# Explicit backend-tier gate. Ordinary package checks exercise all analytic
# behavior; only tests that require a compiled Stan backend opt into this gate.
skip_unless_bayesim_backend <- function() {
  tier <- tolower(Sys.getenv("BAYESIM_TEST_TIER", "core"))
  testthat::skip_if_not(
    tier %in% c("backend", "full"),
    "requires the explicit bayesim backend test tier"
  )
}

# mirai daemon lifecycle ----------------------------------------------------

# Wait for mirai daemon processes to finish exiting after mirai::daemons(0).
#
# daemons(0) signals the daemons and sleeps a fixed 200 ms; the daemon R
# processes then shut down asynchronously, running end-of-session finalizers.
# Under covr every process that loaded the instrumented package dumps a
# coverage trace on exit (covr:::save_trace), and package_coverage() reads
# each trace file with readRDS as soon as testthat finishes. A daemon still
# writing its trace at that point yields a truncated read ("error reading
# from connection" in covr's merge_coverage) and fails CI — the race is
# tightest in test-transport-purrr-mirai.R, the alphabetically last test file,
# whose teardown fires seconds before the traces are read.
#
# Polls ps for live (non-zombie) processes whose command line matches
# mirai::daemon(<url>). Best effort: no-ops on non-unix (Windows never runs
# under covr here) and when ps is unavailable; bounds the wait at timeout.
wait_for_mirai_daemons_exit <- function(timeout = 30, poll = 0.05) {
  if (.Platform$OS.type != "unix" || Sys.which("ps") == "") {
    return(invisible(FALSE))
  }
  ps_lines <- function() {
    tryCatch(
      suppressWarnings(system("ps -eo pid=,stat=,args=", intern = TRUE)),
      error = function(e) character()
    )
  }
  n_live_daemons <- function() {
    # Zombies (stat Z) have exited and can no longer write trace files.
    length(Filter(
      function(l) {
        grepl("mirai::daemon(", l, fixed = TRUE) &&
          !grepl("^[0-9]+ +Z", l)
      },
      ps_lines()
    ))
  }
  deadline <- Sys.time() + timeout
  repeat {
    if (n_live_daemons() == 0L || Sys.time() >= deadline) {
      break
    }
    Sys.sleep(poll)
  }
  invisible(TRUE)
}

# Set mirai daemons for the duration of the current test and tear them down
# on exit, waiting for the daemon processes to fully exit (see
# wait_for_mirai_daemons_exit for why the wait matters under covr). When
# daemons exit promptly (the common case) the wait is a single ps call.
local_mirai_daemons <- function(n = 2L, timeout = 30) {
  mirai::daemons(n)
  withr::defer(
    {
      mirai::daemons(0)
      wait_for_mirai_daemons_exit(timeout = timeout)
    },
    envir = parent.frame()
  )
  invisible(NULL)
}
