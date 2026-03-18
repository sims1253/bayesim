#' @title Utility Functions for bayesim
#' @description Internal utility functions for atomic file operations, hashing,
#'   checksums, timing, error capture, and task ID formatting.
#' @name utils
#' @keywords internal
NULL

# Operators ----------------------------------------------------------------

#' String concatenation operator
#'
#' @param x First string
#' @param y Second string
#' @return Concatenated string
#'
#' @name string-concat
#' @export
#' @keywords internal
`%+%` <- function(x, y) {
  paste0(x, y)
}

#' Null coalescing operator
#'
#' Returns `x` if not NULL, otherwise `y`.
#'
#' @param x Value to check
#' @param y Default value if x is NULL
#' @return `x` if not NULL, otherwise `y`
#'
#' @name null-coalescing
#' @export
#'
#' @examples
#' NULL %||% "default"  # returns "default"
#' "value" %||% "default"  # returns "value"
`%||%` <- function(x, y) if (is.null(x)) y else x

# Atomic File Operations --------------------------------------------------

#' Write JSON file atomically
#'
#' Writes a JSON file using an atomic write operation to prevent corruption
#' from interrupted writes. Writes to a temporary file first, then renames.
#'
#' @param x Object to write to JSON
#' @param path Character string giving the file path
#'
#' @return Invisible NULL. Called for side effect.
#'
#' @details
#' The function writes to a `.tmp` suffix file first, then uses `file.rename()`
#' to atomically move it to the target path. If the rename fails, a
#' `bayesim_checkpoint_error` is thrown.
#'
#' @keywords internal
#' @export
write_json_atomic <- function(x, path) {
  tmp_path <- paste0(path, ".tmp")

  jsonlite::write_json(x, tmp_path, auto_unbox = TRUE, pretty = TRUE)

  success <- file.rename(tmp_path, path)
  if (!success) {
    stop(bayesim_checkpoint_error(
      paste0("Failed to atomically write JSON file to: ", path)
    ))
  }

  invisible(NULL)
}

#' Write RDS file atomically
#'
#' Writes an RDS file using an atomic write operation to prevent corruption
#' from interrupted writes. Writes to a temporary file first, then renames.
#'
#' @param x Object to write to RDS
#' @param path Character string giving the file path
#'
#' @return Invisible NULL. Called for side effect.
#'
#' @keywords internal
#' @export
write_rds_atomic <- function(x, path) {
  tmp_path <- paste0(path, ".tmp")

  saveRDS(x, tmp_path)

  success <- file.rename(tmp_path, path)
  if (!success) {
    stop(bayesim_checkpoint_error(
      paste0("Failed to atomically write RDS file to: ", path)
    ))
  }

  invisible(NULL)
}

# Hashing -----------------------------------------------------------------

#' Compute hash of an R object
#'
#' Computes a hash of an R object using xxHash64 algorithm for fast,
#' deterministic hashing.
#'
#' @param x Any R object to hash
#'
#' @return Character string containing the hash value.
#'
#' @keywords internal
#' @export
compute_hash <- function(x) {
  digest::digest(x, algo = "xxhash64")
}

# Checksums ---------------------------------------------------------------

#' Compute checksum for a file
#'
#' Computes an MD5 checksum for a file. Useful for verifying file integrity.
#'
#' @param path Character string giving the file path
#'
#' @return Character string containing the MD5 checksum.
#'
#' @keywords internal
#' @export
compute_file_checksum <- function(path) {
  digest::digest(file = path, algo = "md5")
}

#' Write checksums for multiple files
#'
#' Computes checksums for specified files in a directory and writes them
#' to a checksums.json file.
#'
#' @param dir_path Character string giving the directory path
#' @param files Character vector of file names (relative to dir_path)
#'
#' @return Invisible NULL. Called for side effect.
#'
#' @keywords internal
#' @export
write_checksums <- function(dir_path, files) {
  checksums <- list()

  for (file in files) {
    file_path <- file.path(dir_path, file)
    if (file.exists(file_path)) {
      checksums[[file]] <- compute_file_checksum(file_path)
    }
  }

  write_json_atomic(checksums, file.path(dir_path, "checksums.json"))

  invisible(NULL)
}

#' Verify checksums for files in a directory
#'
#' Reads the checksums.json file and verifies all files match their recorded
#' checksums.
#'
#' @param dir_path Character string giving the directory path
#'
#' @return Logical. TRUE if all checksums match, FALSE otherwise.
#'
#' @keywords internal
#' @export
verify_checksums <- function(dir_path) {
  checksums_path <- file.path(dir_path, "checksums.json")

  if (!file.exists(checksums_path)) {
    return(FALSE)
  }

  stored <- jsonlite::read_json(checksums_path)

  for (file in names(stored)) {
    file_path <- file.path(dir_path, file)

    if (!file.exists(file_path)) {
      return(FALSE)
    }

    current_checksum <- compute_file_checksum(file_path)
    if (current_checksum != stored[[file]]) {
      return(FALSE)
    }
  }

  TRUE
}

# Timing ------------------------------------------------------------------

#' Create a timer object
#'
#' Creates a timer object with methods for tracking elapsed time across
#' execution phases.
#'
#' @return A list with the following methods:
#'   \itemize{
#'     \item `start()` - Start or restart the timer
#'     \item `stop()` - Stop the timer
#'     \item `elapsed()` - Get elapsed time in seconds
#'   }
#'
#' @export
#' @keywords internal
#' @examples
#' timer <- make_timer()
#' timer$start()
#' Sys.sleep(0.1)
#' timer$stop()
#' timer$elapsed()
make_timer <- function() {
  start_time <- NULL
  stop_time <- NULL

  list(
    start = function() {
      start_time <<- Sys.time()
      stop_time <<- NULL
      invisible(NULL)
    },

    stop = function() {
      stop_time <<- Sys.time()
      invisible(NULL)
    },

    elapsed = function() {
      if (is.null(start_time)) {
        return(0)
      }

      end <- if (is.null(stop_time)) Sys.time() else stop_time

      as.numeric(difftime(end, start_time, units = "secs"))
    }
  )
}

# Error Capture -----------------------------------------------------------

#' Capture error information
#'
#' Captures detailed information about an error condition, including
#' class, message, call, and a trimmed traceback.
#'
#' @param e A condition object (typically from tryCatch)
#'
#' @return A named list with elements:
#'   \itemize{
#'     \item `error_class` - The class of the error
#'     \item `error_message` - The error message
#'     \item `call` - The call that caused the error (if available)
#'     \item `traceback` - Trimmed traceback (limited to 20 frames)
#'   }
#'
#' @keywords internal
#' @export
capture_error_info <- function(e) {
  # Get traceback and trim it
  tb <- tryCatch(
    sys.calls(),
    error = function(cond) NULL
  )

  # Trim traceback to last 20 frames to avoid excessive storage
  if (!is.null(tb) && length(tb) > 20) {
    tb <- tail(tb, 20)
  }

  # Convert calls to character strings
  tb_str <- if (!is.null(tb)) {
    vapply(tb, function(x) deparse(x)[1], character(1))
  } else {
    character(0)
  }

  # Get call if available
  error_call <- tryCatch(
    {
      cl <- conditionCall(e)
      if (is.null(cl)) NULL else deparse(cl)[1]
    },
    error = function(cond) NULL
  )

  list(
    error_class = paste(class(e), collapse = ", "),
    error_message = conditionMessage(e),
    call = error_call,
    traceback = tb_str
  )
}

# Task ID Formatting ------------------------------------------------------

#' Format task ID from indices (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Use [make_task_id()] instead, which auto-calculates field widths
#' from the grid dimensions.
#'
#' @param data_idx Integer. Data index (1-999)
#' @param fit_idx Integer. Fit index (1-999)
#' @param rep_idx Integer. Replication index (1-99999)
#'
#' @return Character string in format "dXXX_fXXX_rXXXXX"
#'
#' @export
#' @examples
#' # Migration: replace format_task_id() with make_task_id()
#' # Before (deprecated):
#' #   format_task_id(1, 2, 100)
#' # After (preferred):
#' make_task_id(1, 2, 100)
format_task_id <- function(data_idx, fit_idx, rep_idx) {
  lifecycle::deprecate_warn(
    "1.1",
    "format_task_id()",
    "make_task_id()"
  )
  make_task_id(
    data_idx,
    fit_idx,
    rep_idx,
    widths = list(data = 3, fit = 3, rep = 5)
  )
}

# Flatten Nested List -----------------------------------------------------

#' Flatten a nested list with prefix
#'
#' Flattens a named list with nested named numeric vectors, using a
#' prefix__name__subname naming convention for the flattened elements.
#'
#' @param x A named list, possibly containing nested named numeric vectors
#' @param prefix Character string to use as prefix for flattened names
#'
#' @return A flattened named list where nested named numeric vectors are
#'   expanded. When `prefix` is non-empty, flattened names use
#'   "prefix__name__sub_name"; when `prefix` is empty (""), they use
#'   "name__sub_name". Scalar and unnamed elements are passed through unchanged.
#'
#' @details
#' This function handles lists where some elements may be named numeric
#' vectors. Those vectors are flattened into individual scalar elements
#' with a double-underscore naming convention. The empty-prefix variant
#' is used internally by checkpointing code where the metric name already
#' serves as the outer namespace.
#'
#' @export
#' @keywords internal
#' @examples
#' x <- list(a = 1, b = c(x = 2, y = 3), c = 4)
#' flatten_with_prefix(x, "param")
#' # Returns: list(a = 1, param__b__x = 2, param__b__y = 3, c = 4)
flatten_with_prefix <- function(x, prefix) {
  result <- list()

  for (name in names(x)) {
    value <- x[[name]]

    if (is.numeric(value) && length(value) > 1 && !is.null(names(value))) {
      # Flatten named numeric vector
      for (sub_name in names(value)) {
        new_name <- if (prefix == "") {
          paste0(name, "__", sub_name)
        } else {
          paste0(prefix, "__", name, "__", sub_name)
        }
        result[[new_name]] <- value[[sub_name]]
      }
    } else if (prefix != "" && length(value) == 1) {
      # Add prefix for scalars when a non-empty prefix is provided
      result[[paste0(prefix, "__", name)]] <- value
    } else {
      result[[name]] <- value
    }
  }

  result
}
