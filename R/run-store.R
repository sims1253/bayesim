#' Internal run-store seam
#'
#' The run store keeps the latest checkpoint state between writes and uses
#' the checkpoint functions to read and persist outcomes. Runs without a
#' result path keep outcomes in the execution loop and do not use a store.
#' @name run-store-internals
#' @keywords internal
NULL

#' Create a filesystem run store.
#'
#' @param result_path Filesystem directory, or NULL to disable checkpoint storage.
#' @param config_fingerprint Study fingerprint.
#' @param config_spec Optional manifest specification.
#' @param checkpoint_format Checkpoint serialization format.
#' @param keep_checkpoints Number of checkpoint commit directories to retain.
#'   Pruning removes old commit directories only; immutable outcome shards and
#'   ledger history are never pruned, so durable storage grows roughly
#'   linearly with completed tasks.
#' @return An internal run-store object with initialize/read/write methods,
#'   or NULL when result_path is NULL.
#' @keywords internal
new_run_store <- function(
  result_path = NULL,
  config_fingerprint = NULL,
  config_spec = NULL,
  checkpoint_format = "rds",
  keep_checkpoints = 2L,
  retention_spec = NULL,
  run_policy_spec = NULL
) {
  if (is.null(result_path)) {
    return(NULL)
  }

  state <- new.env(parent = emptyenv())
  state$checkpoint <- NULL
  store <- list(
    backend = "filesystem",
    path = result_path,
    initialize = function() {
      init_checkpoint_dir(
        result_path,
        config_fingerprint = config_fingerprint,
        config_spec = config_spec,
        checkpoint_format = checkpoint_format,
        retention_spec = retention_spec,
        run_policy_spec = run_policy_spec
      )
      invisible(TRUE)
    },
    read = function() {
      state$checkpoint <- get_latest_valid_checkpoint(
        result_path,
        config_fingerprint = config_fingerprint
      )
      state$checkpoint
    },
    write = function(
      task_grid,
      task_results,
      prior_results_df = NULL,
      prior_task_results = NULL,
      adaptive_next_check = NULL,
      adaptive_state = NULL
    ) {
      state$checkpoint <- write_checkpoint(
        result_path,
        task_grid,
        task_results,
        config_fingerprint = config_fingerprint,
        checkpoint_format = checkpoint_format,
        keep_checkpoints = keep_checkpoints,
        prior_results_df = prior_results_df,
        prior_task_results = prior_task_results,
        adaptive_next_check = adaptive_next_check,
        adaptive_state = adaptive_state,
        run_policy_spec = run_policy_spec,
        prior_checkpoint = state$checkpoint,
        return_state = TRUE,
        delta_store = TRUE
      )
      invisible(state$checkpoint$checkpoint_id)
    }
  )
  class(store) <- "bayesim_run_store"
  store
}

is_run_store <- function(x) inherits(x, "bayesim_run_store")
