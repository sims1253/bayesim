#' @title Bayesim Error Conditions
#' @description Error condition classes for the bayesim package.
#'
#' All errors inherit from `bayesim_error`. Fatal errors (config, contract,
#' checkpoint, internal, validation) stop the simulation run. Recoverable errors
#' (data, fit, metric) allow remaining tasks to continue.
#'
#' @name bayesim-errors
#' @keywords internal
NULL

#' Base Bayesim Error
#'
#' @param message The error message.
#' @param call The call that caused the error (optional).
#' @return An error condition object.
#' @export
bayesim_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_error", "error", "condition")
  )
}

#' Check if a condition is a bayesim error
#'
#' @param cond A condition object to test.
#' @return TRUE if the condition is a bayesim error, FALSE otherwise.
#' @export
is_bayesim_error <- function(cond) {
  inherits(cond, "bayesim_error")
}

# ============================================================================
# FATAL ERRORS
# ============================================================================

#' Configuration validation error (fatal)
#'
#' @inheritParams bayesim_error
#' @export
bayesim_config_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_config_error", "bayesim_error", "error", "condition")
  )
}

#' Contract/interface violation error (fatal)
#'
#' @inheritParams bayesim_error
#' @export
bayesim_contract_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_contract_error", "bayesim_error", "error", "condition")
  )
}

#' Checkpoint read/write error (fatal)
#'
#' @inheritParams bayesim_error
#' @export
bayesim_checkpoint_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_checkpoint_error", "bayesim_error", "error", "condition")
  )
}

#' Internal consistency error (fatal)
#'
#' @inheritParams bayesim_error
#' @export
bayesim_internal_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_internal_error", "bayesim_error", "error", "condition")
  )
}

#' Validation error (fatal)
#'
#' Inherits from `bayesim_contract_error`.
#'
#' @inheritParams bayesim_error
#' @export
bayesim_validation_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c(
      "bayesim_validation_error",
      "bayesim_contract_error",
      "bayesim_error",
      "error",
      "condition"
    )
  )
}

# ============================================================================
# RECOVERABLE ERRORS
# ============================================================================

#' Data generation/validation error (recoverable)
#'
#' @inheritParams bayesim_error
#' @export
bayesim_data_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_data_error", "bayesim_error", "error", "condition")
  )
}

#' Model fitting error (recoverable)
#'
#' @inheritParams bayesim_error
#' @export
bayesim_fit_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_fit_error", "bayesim_error", "error", "condition")
  )
}

#' Metric computation error (recoverable)
#'
#' @inheritParams bayesim_error
#' @export
bayesim_metric_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_metric_error", "bayesim_error", "error", "condition")
  )
}

# ============================================================================
# ERROR CLASSIFICATION HELPERS
# ============================================================================

#' Check if an error is fatal
#'
#' @param cond A condition object to test.
#' @return TRUE if the error should stop the entire simulation run.
#' @export
is_fatal_error <- function(cond) {
  inherits(cond, "bayesim_config_error") ||
    inherits(cond, "bayesim_contract_error") ||
    inherits(cond, "bayesim_checkpoint_error") ||
    inherits(cond, "bayesim_internal_error") ||
    inherits(cond, "bayesim_validation_error")
}

#' Check if an error is recoverable (task-level)
#'
#' @param cond A condition object to test.
#' @return TRUE if the simulation can continue with other tasks.
#' @export
is_recoverable_error <- function(cond) {
  inherits(cond, "bayesim_data_error") ||
    inherits(cond, "bayesim_fit_error") ||
    inherits(cond, "bayesim_metric_error")
}
