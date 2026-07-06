#' @title S7 Property Validator Helpers
#' @description Factory functions that return validator closures suitable for
#'   use as the `validator =` argument of [S7::new_property()]. Each factory
#'   captures a small, commonly-repeated validation rule and returns a
#'   `function(value)` that follows the S7 convention: it returns `NULL` when
#'   `value` is valid, or a single character string describing why `value` is
#'   invalid.
#'
#'   These helpers exist to deduplicate the small S7 property validators that
#'   are otherwise copy-pasted across class definitions (B5).
#' @name s7-validators
#' @keywords internal
NULL

# =============================================================================
# Validator factories
# =============================================================================

#' Validate a Positive Integer
#'
#' Returns an S7 property validator function that requires `value` to be a
#' single non-NA integer greater than or equal to 1. The returned function
#' returns `NULL` when valid, otherwise the string `"<message>"` (defaults to
#' `"must be a positive integer"`).
#'
#' @param message Character scalar; the error string returned when `value` is
#'   invalid. Defaults to `"must be a positive integer"`. Pass a
#'   property-name-prefixed string to preserve a previously established error
#'   message verbatim.
#' @return A `function(value)` suitable as an S7 property `validator`.
#' @keywords internal
validate_positive_integer <- function(message = "must be a positive integer") {
  force(message)
  function(value) {
    if (length(value) != 1L || is.na(value) || value < 1) {
      message
    }
  }
}
