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

#' Validate a Single Character String
#'
#' Returns an S7 property validator function that requires `value` to be a
#' length-1 character vector. When `allow_na = FALSE` (the default) an NA
#' value is rejected. The returned function returns `NULL` when valid,
#' otherwise the string `"must be a single non-NA character string"` (or, when
#' `allow_na = TRUE`, `"must be a single character string"`).
#'
#' @param allow_na Logical; if `TRUE`, NA values are accepted.
#' @return A `function(value)` suitable as an S7 property `validator`.
#' @keywords internal
validate_scalar_character <- function(allow_na = FALSE) {
  function(value) {
    if (length(value) != 1L || (is.na(value) && !allow_na)) {
      if (allow_na) {
        "must be a single character string"
      } else {
        "must be a single non-NA character string"
      }
    }
  }
}

#' Validate a Positive Numeric
#'
#' Returns an S7 property validator function that requires `value` to be a
#' single non-NA numeric greater than or equal to zero. The returned function
#' returns `NULL` when valid, otherwise the string
#' `"must be a positive numeric"`.
#'
#' @return A `function(value)` suitable as an S7 property `validator`.
#' @keywords internal
validate_positive_numeric <- function() {
  function(value) {
    if (length(value) != 1L || is.na(value) || value < 0) {
      "must be a positive numeric"
    }
  }
}

#' Validate a Character Vector With No NA / Empty Strings
#'
#' Returns an S7 property validator function that requires `value` to be a
#' character vector containing no `NA` and no empty (`""`) elements. The
#' returned function returns `NULL` when valid, otherwise the string
#' `"must not contain NA or empty strings"`.
#'
#' @return A `function(value)` suitable as an S7 property `validator`.
#' @keywords internal
validate_character_no_na <- function() {
  function(value) {
    if (anyNA(value) || any(value == "")) {
      "must not contain NA or empty strings"
    }
  }
}

#' Validate a Vector Is a Subset of Allowed Values
#'
#' Returns an S7 property validator function that requires every element of
#' `value` to be present in `allowed`. The returned function returns `NULL`
#' when valid, otherwise a string of the form
#' `"must be a subset of: <allowed>"`.
#'
#' @param allowed A character vector of permitted values.
#' @return A `function(value)` suitable as an S7 property `validator`.
#' @keywords internal
validate_subset <- function(allowed) {
  function(value) {
    bad <- setdiff(value, allowed)
    if (length(bad) > 0L) {
      paste0(
        "must be a subset of: ",
        paste(allowed, collapse = ", ")
      )
    }
  }
}
