#' @import stats
#' @importFrom utils capture.output tail globalVariables
#' @importFrom withr with_seed
#' @importFrom loo loo psis E_loo relative_eff pareto_k_values
#' @importFrom posterior ndraws as_draws_matrix as_draws_df ess_bulk
#' @keywords internal
utils::globalVariables("mb")
NULL
