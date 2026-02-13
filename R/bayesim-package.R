#' @keywords internal
#' @import rlang
#' @importFrom cli cli_abort cli_warn cli_inform
#' @importFrom lifecycle deprecated
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' Simulation Framework for Bayesian Modeling
#'
#' @description
#' bayesim provides a modern simulation framework for Bayesian modeling studies.
#' It offers extensible tools for running complex simulation studies with
#' deterministic reproducibility, checkpoint/resume capabilities, and
#' memory-bounded execution.
#'
#' @details
#' Key features:
#' \itemize{
#'   \item Extensible S7-based interfaces for custom fitters and metrics
#'   \item Deterministic reproducibility across sequential/parallel/resume modes
#'   \item File-based checkpoint and resume with atomic writes
#'   \item Memory-bounded execution with configurable artifact retention
#' }
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item [simulation_config()] - Configure simulation runs
#'   \item [run_simulation()] - Execute simulation
#'   \item [Fitter] - Interface for model fitters
#'   \item [Metric] - Interface for metrics
#' }
#'
#' @docType package
#' @name bayesim-package
NULL
