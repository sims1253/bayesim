#' Internal lifecycle contracts
#'
#' The public simulation API intentionally keeps the historical `skipped`
#' status for compatibility.  A skipped task is not, however, a completed
#' scientific outcome: `stop_reason` records why execution stopped and resume
#' turns policy-stopped tasks back into pending work.
#' @name lifecycle-internals
#' @keywords internal
NULL

TASK_STATUSES <- c(
  "pending",
  "success",
  "failed",
  "skipped"
)

TASK_TERMINAL_STATUSES <- c("success", "failed")
TASK_RESUMABLE_STATUSES <- setdiff(TASK_STATUSES, TASK_TERMINAL_STATUSES)

is_valid_task_status <- function(status) {
  if (!is.character(status)) {
    return(FALSE)
  }
  !is.na(status) & status %in% TASK_STATUSES
}

#' Return whether a task status represents a completed scientific outcome.
#' @keywords internal
is_terminal_task_status <- function(status) {
  if (!is.character(status)) {
    return(FALSE)
  }
  !is.na(status) & status %in% TASK_TERMINAL_STATUSES
}

#' Return whether a task status is eligible for execution on resume.
#' @keywords internal
is_resumable_task_status <- function(status) {
  if (!is.character(status)) {
    return(FALSE)
  }
  !is.na(status) & status %in% TASK_RESUMABLE_STATUSES
}

#' Return whether a result path has no user data yet.
#' @keywords internal
is_empty_result_path <- function(result_path) {
  if (!dir.exists(result_path)) {
    return(!file.exists(result_path))
  }
  length(list.files(result_path, all.files = TRUE, no.. = TRUE)) == 0L
}

#' Classify a result path before a run starts.
#' @keywords internal
result_path_state <- function(result_path) {
  if (is.null(result_path)) {
    return("memory")
  }
  if (file.exists(result_path) && !dir.exists(result_path)) {
    return("file")
  }
  if (!dir.exists(result_path)) {
    return("missing")
  }
  if (is_empty_result_path(result_path)) {
    return("empty")
  }
  if (file.exists(file.path(result_path, "run_manifest.json"))) {
    return("run")
  }
  "ambiguous"
}

#' Construct the immutable scientific part of a simulation configuration.
#' @keywords internal
new_study_spec <- function(config) {
  if (!is_simulation_config(config)) {
    cli::cli_abort("config must be a SimulationConfig object")
  }
  structure(
    list(
      data_grid = config@data_grid,
      fit_grid = config@fit_grid,
      task_grid = config@task_grid,
      data_generator = config@data_generator,
      fitter = config@fitter,
      metrics = config@metrics,
      n_replicates = config@n_replicates,
      seed = config@seed
    ),
    class = c("bayesim_study_spec", "list")
  )
}

#' Construct the operational policy for a simulation run.
#' @keywords internal
new_run_policy <- function(
  config,
  progress = TRUE,
  workers = NULL,
  verbose = TRUE
) {
  if (!is_simulation_config(config)) {
    cli::cli_abort("config must be a SimulationConfig object")
  }
  structure(
    list(
      result_path = config@result_path,
      checkpoint_every = config@checkpoint_every,
      checkpoint_format = config@checkpoint_format,
      keep_checkpoints = config@keep_checkpoints,
      retain = config@retain,
      max_errors = config@max_errors,
      stop_on = config@stop_on,
      progress = progress,
      verbose = verbose,
      workers = workers
    ),
    class = c("bayesim_run_policy", "list")
  )
}

is_study_spec <- function(x) inherits(x, "bayesim_study_spec")
is_run_policy <- function(x) inherits(x, "bayesim_run_policy")

#' Convert operational policy to a manifest-safe resume specification.
#' @keywords internal
as_run_policy_spec <- function(policy) {
  if (!is_run_policy(policy)) {
    cli::cli_abort("policy must be a bayesim RunPolicy")
  }
  list(
    checkpoint_every = policy$checkpoint_every,
    checkpoint_format = policy$checkpoint_format,
    keep_checkpoints = policy$keep_checkpoints,
    retain = policy$retain,
    max_errors = policy$max_errors,
    stop_on = policy$stop_on
  )
}
