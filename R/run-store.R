#' Internal run-store seam
#'
#' The execution engine currently uses the established checkpoint functions,
#' but all durable state can be addressed through this small adapter. Keeping
#' this seam separate makes lifecycle tests independent of the filesystem and
#' gives future append-only stores a single contract to implement.
#' @name run-store-internals
#' @keywords internal
NULL

#' Create a run store backed by memory or the filesystem.
#'
#' @param result_path NULL for memory, otherwise a filesystem directory.
#' @param config_fingerprint Study fingerprint.
#' @param config_spec Optional manifest specification.
#' @param checkpoint_format Checkpoint serialization format.
#' @param keep_checkpoints Number of checkpoint commit directories to retain.
#'   Pruning removes old commit directories only; immutable outcome shards and
#'   ledger history are never pruned, so durable storage grows roughly
#'   linearly with completed tasks.
#' @return An internal run-store object with initialize/read/write methods.
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
    state <- new.env(parent = emptyenv())
    state$checkpoint_id <- 0L
    state$checkpoint <- NULL
    store <- list(
      backend = "memory",
      path = NULL,
      initialize = function() invisible(TRUE),
      read = function() state$checkpoint,
      write = function(
        task_grid,
        task_results,
        prior_results_df = NULL,
        prior_task_results = NULL,
        adaptive_next_check = NULL,
        adaptive_state = NULL
      ) {
        state$checkpoint_id <- state$checkpoint_id + 1L
        current_outcomes <- Filter(
          function(x) !is.null(x) && is_bayesim_task_result(x),
          task_results
        )
        prior_outcomes <- Filter(
          function(x) !is.null(x) && is_bayesim_task_result(x),
          prior_task_results %||% list()
        )
        outcomes <- prior_outcomes
        outcome_ids <- if (length(outcomes) > 0L) {
          vapply(outcomes, function(x) x$task_id, character(1))
        } else {
          character()
        }
        # Assigning current outcomes last makes a store write an atomic
        # replacement by task identity, matching the filesystem adapter.
        for (outcome in current_outcomes) {
          hit <- match(outcome$task_id, outcome_ids)
          if (is.na(hit)) {
            outcomes <- c(outcomes, list(outcome))
            outcome_ids <- c(outcome_ids, outcome$task_id)
          } else {
            outcomes[[hit]] <- outcome
          }
        }
        state$checkpoint <- list(
          checkpoint_id = state$checkpoint_id,
          task_grid = task_grid,
          task_outcomes = outcomes,
          results_df = results_to_dataframe(outcomes),
          adaptive_next_check = adaptive_next_check,
          adaptive_state = adaptive_state,
          meta = list(
            adaptive_next_check = adaptive_next_check,
            adaptive_state = adaptive_state,
            run_policy = run_policy_spec
          )
        )
        invisible(state$checkpoint_id)
      }
    )
    class(store) <- "bayesim_run_store"
    return(store)
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
