#' bayesim: Simulation Framework for Bayesian Modeling
#'
#' bayesim provides a simulation framework for Bayesian modeling studies with
#' reproducible execution, checkpoint/resume support, and memory-bounded task
#' processing.
#'
#' @keywords internal
#' @import rlang
#' @importFrom cli cli_abort cli_warn cli_inform
#' @importFrom lifecycle deprecated deprecate_warn
#' @seealso [simulation_config()], [run_simulation()], [resume_simulation()],
#'   [Fitter], [Metric]
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

# ============================================================================
# Schema Version Constants
# ============================================================================

#' Run Schema Version
#'
#' Version identifier for checkpoint format compatibility.
#' Increment this when the on-disk checkpoint format changes in a way
#' that breaks backward compatibility.
#'
#' @keywords internal
RUN_SCHEMA_VERSION <- 1L

#' Result Schema Version
#'
#' Version identifier for result column contract.
#' Increment this when result column names or types change.
#'
#' @keywords internal
RESULT_SCHEMA_VERSION <- 1L
