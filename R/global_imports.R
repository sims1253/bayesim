#' @import stats
#' @importFrom utils capture.output tail globalVariables
#' @importFrom withr with_seed
#' @importFrom loo loo psis E_loo relative_eff pareto_k_values
#' @importFrom posterior ndraws as_draws_matrix as_draws_df ess_bulk
#' @keywords internal
utils::globalVariables("mb")
NULL

# carrier is used at runtime by purrr::in_parallel() (it crates the task
# lambda shipped to mirai daemons); this reference keeps the declared Import
# visible to R CMD check.
# jarl-ignore unused_function: deliberate no-op reference for R CMD check
.reference_unused_imports <- function() {
  carrier::crate
}
