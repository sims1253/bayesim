#' @title Bayesim Error Conditions
#' @description Error condition classes for the bayesim package. These provide
#'   structured error handling with explicit classification of error severity
#'   and type.
#' @name bayesim-errors
#' @keywords internal
NULL

#' @title Base Bayesim Error
#' @description Base class for all bayesim errors. All bayesim-specific errors
#'   inherit from this class in addition to "error" and "condition".
#' @param message The error message
#' @param call The call that caused the error (optional)
#' @return An error condition object
#' @keywords internal
#' @export
bayesim_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_error", "error", "condition")
  )
}

#' Check if a condition is a bayesim error
#'
#' @param cond A condition object to test
#' @return TRUE if the condition is a bayesim error, FALSE otherwise
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' tryCatch(
#'   bayesim_config_error("Invalid config"),
#'   bayesim_error = function(cond) {
#'     is_bayesim_error(cond) # TRUE
#'   }
#' )
#' }
is_bayesim_error <- function(cond) {
  inherits(cond, "bayesim_error")
}

# ============================================================================
# FATAL ERRORS - Stop entire simulation run
# ============================================================================

#' Configuration validation error (Fatal)
#'
#' @description
#' Raised when configuration parameters are invalid, missing, or inconsistent.
#' This is a fatal error that stops the entire simulation run.
#'
#' @inheritParams bayesim_error
#' @return An error condition object with class c("bayesim_config_error",
#'   "bayesim_error", "error", "condition")
#' @keywords internal
#' @export
#' @seealso [is_fatal_error()]
bayesim_config_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_config_error", "bayesim_error", "error", "condition")
  )
}

#' Contract/interface violation error (Fatal)
#'
#' @description
#' Raised when a function contract or interface is violated, such as
#' incorrect argument types or unexpected return values. This is a fatal
#' error that stops the entire simulation run.
#'
#' @inheritParams bayesim_error
#' @return An error condition object with class c("bayesim_contract_error",
#'   "bayesim_error", "error", "condition")
#' @keywords internal
#' @export
#' @seealso [is_fatal_error()]
bayesim_contract_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_contract_error", "bayesim_error", "error", "condition")
  )
}

#' Checkpoint read/write error (Fatal)
#'
#' @description
#' Raised when checkpoint file operations fail, such as reading corrupted
#' checkpoint data or failing to write checkpoint files. This is a fatal
#' error that stops the entire simulation run.
#'
#' @inheritParams bayesim_error
#' @return An error condition object with class c("bayesim_checkpoint_error",
#'   "bayesim_error", "error", "condition")
#' @keywords internal
#' @export
#' @seealso [is_fatal_error()]
bayesim_checkpoint_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_checkpoint_error", "bayesim_error", "error", "condition")
  )
}

#' Internal consistency error (Fatal)
#'
#' @description
#' Raised when an internal consistency check fails. These errors indicate
#' bugs in the bayesim code itself and should never occur in correct code.
#' This is a fatal error that stops the entire simulation run.
#'
#' @inheritParams bayesim_error
#' @return An error condition object with class c("bayesim_internal_error",
#'   "bayesim_error", "error", "condition")
#' @keywords internal
#' @export
#' @seealso [is_fatal_error()]
bayesim_internal_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_internal_error", "bayesim_error", "error", "condition")
  )
}

#' Validation error (Fatal)
#'
#' @description
#' Raised when input validation fails. This is a fatal error that stops
#' the entire simulation run.
#'
#' @inheritParams bayesim_error
#' @return An error condition object with class c("bayesim_validation_error",
#'   "bayesim_contract_error", "bayesim_error", "error", "condition")
#' @keywords internal
#' @export
#' @seealso [is_fatal_error()]
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
# RECOVERABLE ERRORS - Task-level, can continue with other tasks
# ============================================================================

#' Data generation/validation error (Recoverable)
#'
#' @description
#' Raised when data generation or validation fails for a specific task.
#' This is a task-level recoverable error that allows the simulation to
#' continue with other tasks.
#'
#' @inheritParams bayesim_error
#' @return An error condition object with class c("bayesim_data_error",
#'   "bayesim_error", "error", "condition")
#' @keywords internal
#' @export
#' @seealso [is_recoverable_error()]
bayesim_data_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_data_error", "bayesim_error", "error", "condition")
  )
}

#' Model fitting error (Recoverable)
#'
#' @description
#' Raised when model fitting fails for a specific task, such as convergence
#' failures or numerical issues. This is a task-level recoverable error that
#' allows the simulation to continue with other tasks.
#'
#' @inheritParams bayesim_error
#' @return An error condition object with class c("bayesim_fit_error",
#'   "bayesim_error", "error", "condition")
#' @keywords internal
#' @export
#' @seealso [is_recoverable_error()]
bayesim_fit_error <- function(message, call = NULL) {
  structure(
    list(message = message, call = call),
    class = c("bayesim_fit_error", "bayesim_error", "error", "condition")
  )
}

#' Metric computation error (Recoverable)
#'
#' @description
#' Raised when metric computation fails for a specific task, such as
#' undefined metrics for edge cases. This is a task-level recoverable error
#' that allows the simulation to continue with other tasks.
#'
#' @inheritParams bayesim_error
#' @return An error condition object with class c("bayesim_metric_error",
#'   "bayesim_error", "error", "condition")
#' @keywords internal
#' @export
#' @seealso [is_recoverable_error()]
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
#' @description
#' Checks whether an error condition is fatal (should stop the entire
#' simulation run). Fatal errors include configuration errors, contract
#' violations, checkpoint failures, and internal consistency errors.
#'
#' @param cond A condition object to test
#' @return TRUE if the error is fatal, FALSE otherwise
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' tryCatch(
#'   bayesim_config_error("Invalid parameter"),
#'   bayesim_error = function(cond) {
#'     if (is_fatal_error(cond)) {
#'       # Stop everything
#'     }
#'   }
#' )
#' }
is_fatal_error <- function(cond) {
  inherits(cond, "bayesim_config_error") ||
    inherits(cond, "bayesim_contract_error") ||
    inherits(cond, "bayesim_checkpoint_error") ||
    inherits(cond, "bayesim_internal_error") ||
    # bayesim_validation_error inherits from bayesim_contract_error,
    # but is listed explicitly for clarity
    inherits(cond, "bayesim_validation_error")
}

#' Check if an error is recoverable (task-level)
#'
#' @description
#' Checks whether an error condition is recoverable at the task level.
#' Recoverable errors include data generation errors, model fitting errors,
#' and metric computation errors. The simulation can continue with other
#' tasks when these errors occur.
#'
#' @param cond A condition object to test
#' @return TRUE if the error is recoverable, FALSE otherwise
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' tryCatch(
#'   bayesim_fit_error("Model did not converge"),
#'   bayesim_error = function(cond) {
#'     if (is_recoverable_error(cond)) {
#'       # Log error and continue with next task
#'     }
#'   }
#' )
#' }
is_recoverable_error <- function(cond) {
  inherits(cond, "bayesim_data_error") ||
    inherits(cond, "bayesim_fit_error") ||
    inherits(cond, "bayesim_metric_error")
}
