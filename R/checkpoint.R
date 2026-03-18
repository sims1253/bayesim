#' @title Checkpoint Functionality for bayesim
#' @description Functions for creating, reading, and managing checkpoints during
#'   simulation runs. Checkpoints provide atomic write operations, integrity
#'   verification via checksums, and support for resume functionality.
#' @name checkpoint
#' @keywords internal
NULL

# =============================================================================
# Schema Version Constants
# =============================================================================

#' Checkpoint schema version constants
#'
#' These constants define the schema versions for checkpoint files. They are
#' incremented when incompatible format changes are made.
#'
#' @name checkpoint-constants
#' @keywords internal
NULL

#' Run Schema Version
#'
#' Version identifier for checkpoint format compatibility.
#' Increment this when the on-disk checkpoint format changes in a way
#' that breaks backward compatibility.
#'
#' @keywords internal
RUN_SCHEMA_VERSION <- 1L

#' Result Schema Version
#'
#' Version identifier for result column contract.
#' Increment this when result column names or types change.
#'
#' @keywords internal
RESULT_SCHEMA_VERSION <- 1L

# =============================================================================
# Checkpoint Directory Initialization
# =============================================================================

#' Initialize checkpoint directory
#'
#' Creates the checkpoint directory structure for a new simulation run.
#' This includes the base result path, checkpoints subdirectory, run manifest,
#' and initial latest pointer.
#'
#' @param result_path Character string giving the base path for results.
#'   If NULL, the function returns NULL immediately (no checkpointing).
#' @param config_fingerprint Character string containing a hash of the
#'   configuration. This is stored in the manifest for validation during
#'   resume operations.
#' @param config_spec Optional normalized configuration spec to persist in the
#'   run manifest for future rehydration.
#' @param checkpoint_format Character scalar naming the checkpoint storage format.
#'
#' @return Invisible path to the result directory, or NULL if result_path is NULL.
#'
#' @details
#' The directory structure created is:
#' \preformatted{
#' result_path/
#' +-- run_manifest.json    # run-level metadata and schema versions
#' +-- latest.json          # pointer to latest valid checkpoint ID
#' +-- checkpoints/         # directory for checkpoint snapshots
#' }
#'
#' The run manifest contains:
#' \itemize{
#'   \item `run_schema_version` - Schema version for the run structure
#'   \item `result_schema_version` - Schema version for result files
#'   \item `config_fingerprint` - Hash of the configuration
#'   \item `created` - Timestamp when the run was created
#' }
#'
#' @seealso [write_checkpoint()], [read_checkpoint()]
#'
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' result_path <- init_checkpoint_dir(
#'   result_path = "/path/to/results",
#'   config_fingerprint = "abc123hash"
#' )
#' }
init_checkpoint_dir <- function(
  result_path,
  config_fingerprint,
  config_spec = NULL,
  checkpoint_format = "rds"
) {
  if (is.null(result_path)) {
    return(NULL)
  }

  assert_supported_checkpoint_format(checkpoint_format)

  # Create base directory
  dir.create(result_path, recursive = TRUE, showWarnings = FALSE)

  # Create checkpoints subdirectory
  dir.create(
    file.path(result_path, "checkpoints"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Write run manifest
  manifest <- list(
    run_schema_version = RUN_SCHEMA_VERSION,
    result_schema_version = RESULT_SCHEMA_VERSION,
    config_fingerprint = config_fingerprint,
    config_spec = config_spec,
    checkpoint_format = checkpoint_format,
    created = as.character(Sys.time())
  )
  write_json_atomic(manifest, file.path(result_path, "run_manifest.json"))

  # Initialize latest pointer (no checkpoints yet)
  write_json_atomic(
    list(checkpoint_id = NULL),
    file.path(result_path, "latest.json")
  )

  invisible(result_path)
}

# =============================================================================
# Checkpoint ID Management
# =============================================================================

#' Get next checkpoint ID
#'
#' Determines the next available checkpoint ID by examining existing
#' checkpoint directories.
#'
#' @param result_path Character string giving the base path for results.
#'
#' @return Integer. The next checkpoint ID (1 if no checkpoints exist).
#'
#' @details
#' Checkpoint directories are named with the format `cp_XXXXXX` where X
#' represents a zero-padded six-digit number. This function scans existing
#' directories and returns the maximum ID + 1.
#'
#' @keywords internal
#' @export
get_next_checkpoint_id <- function(result_path) {
  checkpoints_dir <- file.path(result_path, "checkpoints")

  if (!dir.exists(checkpoints_dir)) {
    return(1L)
  }

  existing <- list.dirs(checkpoints_dir, recursive = FALSE, full.names = FALSE)

  # Only count valid checkpoint directories (cp_XXXXXX format)
  # This ignores .tmp directories from failed write attempts
  valid_dirs <- existing[grepl("^cp_\\d{6}$", existing)]
  if (length(valid_dirs) == 0) {
    return(1L)
  }

  # Extract numeric IDs from directory names (e.g., "cp_000001" -> 1)
  ids <- as.integer(gsub("^cp_", "", valid_dirs))
  max(ids, na.rm = TRUE) + 1L
}

# =============================================================================
# Checkpoint Writing
# =============================================================================

#' Write checkpoint
#'
#' Atomically writes a checkpoint containing the task grid (ledger) and
#' results. Uses a write-temp-then-rename protocol to ensure consistency.
#'
#' @param result_path Character string giving the base path for results.
#'   If NULL, the function returns NULL immediately (no checkpointing).
#' @param task_grid A tibble/data.frame containing the task grid with status
#'   information. Each row represents a task with columns for task_id, status,
#'   and other task metadata.
#' @param task_results A list of `bayesim_task_result` objects containing
#'   results from completed tasks.
#' @param config_fingerprint Character string containing a hash of the
#'   configuration for validation purposes.
#' @param checkpoint_format Character scalar naming the checkpoint storage format.
#'
#' @return Invisible checkpoint ID (integer), or NULL if result_path is NULL.
#'
#' @details
#' The checkpoint directory structure is:
#' \preformatted{
#' checkpoints/
#' +-- cp_000001/
#'     +-- meta.json         # checkpoint metadata
#'     +-- ledger.rds        # task grid with status
#'     +-- results.rds       # metrics + diagnostics per task
#'     +-- checksums.json    # file integrity checksums
#' }
#'
#' Atomic write protocol:
#' \enumerate{
#'   \item Create temporary directory with `.tmp` suffix
#'   \item Write all files to temporary directory
#'   \item Compute and write checksums
#'   \item Validate checksums in temporary directory
#'   \item Rename temporary directory to final name (atomic operation)
#'   \item Update latest.json pointer
#' }
#'
#' If any step fails, the temporary directory is cleaned up and an error
#' is thrown. This ensures that partial writes never appear as valid
#' checkpoints.
#'
#' @seealso [read_checkpoint()], [init_checkpoint_dir()]
#'
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' checkpoint_id <- write_checkpoint(
#'   result_path = "/path/to/results",
#'   task_grid = task_grid,
#'   task_results = task_results,
#'   config_fingerprint = "abc123hash"
#' )
#' }
write_checkpoint <- function(
  result_path,
  task_grid,
  task_results,
  config_fingerprint,
  checkpoint_format = "rds"
) {
  if (is.null(result_path)) {
    return(invisible(NULL))
  }

  assert_supported_checkpoint_format(checkpoint_format)

  # Get next checkpoint ID and create directory names
  checkpoint_id <- get_next_checkpoint_id(result_path)
  checkpoint_name <- sprintf("cp_%06d", checkpoint_id)
  checkpoint_dir <- file.path(result_path, "checkpoints", checkpoint_name)
  tmp_dir <- paste0(checkpoint_dir, ".tmp")

  # Clean up any leftover tmp directory from failed previous attempts
  if (dir.exists(tmp_dir)) {
    unlink(tmp_dir, recursive = TRUE)
  }

  # Create temporary directory
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  # Use tryCatch to ensure cleanup on error
  tryCatch(
    {
      # Write meta.json
      meta <- list(
        checkpoint_id = checkpoint_id,
        created = as.character(Sys.time()),
        config_fingerprint = config_fingerprint,
        n_tasks = nrow(task_grid),
        n_success = sum(task_grid$status == "success", na.rm = TRUE),
        n_failed = sum(task_grid$status == "failed", na.rm = TRUE),
        n_pending = sum(task_grid$status == "pending", na.rm = TRUE)
      )
      write_json_atomic(meta, file.path(tmp_dir, "meta.json"))

      ledger_path <- checkpoint_data_path(tmp_dir, "ledger", checkpoint_format)
      results_path <- checkpoint_data_path(
        tmp_dir,
        "results",
        checkpoint_format
      )

      # Write ledger (task grid)
      write_checkpoint_object(task_grid, ledger_path, checkpoint_format)

      # Convert and write results
      results_df <- results_to_dataframe(task_results)

      write_checkpoint_object(results_df, results_path, checkpoint_format)

      # Read-back validation
      test_grid <- tryCatch(
        read_checkpoint_object(ledger_path, checkpoint_format),
        error = function(e) NULL
      )
      if (is.null(test_grid)) {
        stop(bayesim_checkpoint_error(
          "Checkpoint read-back validation failed for ledger"
        ))
      }

      test_results <- tryCatch(
        read_checkpoint_object(results_path, checkpoint_format),
        error = function(e) NULL
      )
      if (is.null(test_results)) {
        stop(bayesim_checkpoint_error(
          "Checkpoint read-back validation failed for results"
        ))
      }

      # Write checksums for integrity verification
      write_checksums(
        tmp_dir,
        c(
          "meta.json",
          basename(ledger_path),
          basename(results_path)
        )
      )

      # Validate the temporary checkpoint
      if (!verify_checksums(tmp_dir)) {
        stop(bayesim_checkpoint_error(
          "Checkpoint validation failed: checksum mismatch"
        ))
      }

      # Atomic rename of directory
      if (!file.rename(tmp_dir, checkpoint_dir)) {
        stop(bayesim_checkpoint_error(
          paste0(
            "Failed to rename checkpoint directory: ",
            tmp_dir,
            " -> ",
            checkpoint_dir
          )
        ))
      }

      # Update latest pointer
      write_json_atomic(
        list(checkpoint_id = checkpoint_id),
        file.path(result_path, "latest.json")
      )

      invisible(checkpoint_id)
    },
    error = function(e) {
      # Clean up temporary directory on any error
      if (dir.exists(tmp_dir)) {
        unlink(tmp_dir, recursive = TRUE)
      }
      stop(e)
    }
  )
}

# =============================================================================
# Checkpoint Reading
# =============================================================================

#' Read checkpoint
#'
#' Reads a checkpoint by ID, or the latest valid checkpoint if ID is not
#' specified. Verifies checksums before returning data.
#'
#' @param result_path Character string giving the base path for results.
#'   If NULL, the function returns NULL immediately.
#' @param checkpoint_id Integer specifying the checkpoint ID to read.
#'   If NULL (default), reads the latest checkpoint from latest.json.
#'
#' @return A list with the following elements, or NULL if checkpoint not found
#'   or invalid:
#'   \itemize{
#'     \item `checkpoint_id` - Integer ID of the checkpoint
#'     \item `meta` - List of checkpoint metadata
#'     \item `task_grid` - Tibble/data.frame of task grid with status
#'     \item `results_df` - Data frame of task results
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item If checkpoint_id is NULL, read latest.json to get the latest ID
#'   \item Construct checkpoint directory path from the ID
#'   \item Verify directory exists
#'   \item Verify checksums match
#'   \item Read meta.json, ledger.rds, and results.rds
#'   \item Return assembled checkpoint data
#' }
#'
#' If the checkpoint has invalid checksums, a warning is issued and NULL
#' is returned. This allows the caller to fall back to earlier checkpoints.
#'
#' @seealso [write_checkpoint()], [list_checkpoints()], [validate_checkpoint_fingerprint()]
#'
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' # Read latest checkpoint
#' checkpoint <- read_checkpoint("/path/to/results")
#'
#' # Read specific checkpoint
#' checkpoint <- read_checkpoint("/path/to/results", checkpoint_id = 5)
#' }
read_checkpoint <- function(result_path, checkpoint_id = NULL) {
  if (is.null(result_path)) {
    return(NULL)
  }

  # If no checkpoint_id specified, get the latest
  if (is.null(checkpoint_id)) {
    latest_path <- file.path(result_path, "latest.json")

    if (!file.exists(latest_path)) {
      return(NULL)
    }

    latest <- jsonlite::read_json(latest_path)
    checkpoint_id <- latest$checkpoint_id

    if (is.null(checkpoint_id) || length(checkpoint_id) == 0) {
      return(NULL)
    }
  }

  # Construct checkpoint directory path
  checkpoint_name <- sprintf("cp_%06d", checkpoint_id)
  checkpoint_dir <- file.path(result_path, "checkpoints", checkpoint_name)

  if (!dir.exists(checkpoint_dir)) {
    return(NULL)
  }

  # Verify checksums before reading
  if (!verify_checksums(checkpoint_dir)) {
    cli::cli_warn("Checkpoint {checkpoint_id} has invalid checksums")
    return(NULL)
  }

  manifest <- read_run_manifest(result_path)
  checkpoint_format <- manifest$checkpoint_format %||% "rds"
  assert_supported_checkpoint_format(checkpoint_format)

  # Read checkpoint files
  meta <- jsonlite::read_json(file.path(checkpoint_dir, "meta.json"))
  task_grid <- read_checkpoint_object(
    checkpoint_data_path(checkpoint_dir, "ledger", checkpoint_format),
    checkpoint_format
  )
  results_df <- read_checkpoint_object(
    checkpoint_data_path(checkpoint_dir, "results", checkpoint_format),
    checkpoint_format
  )

  list(
    checkpoint_id = checkpoint_id,
    meta = meta,
    task_grid = task_grid,
    results_df = results_df
  )
}

# =============================================================================
# Checkpoint Validation
# =============================================================================

#' Validate checkpoint fingerprint
#'
#' Checks whether a checkpoint's configuration fingerprint matches the
#' expected fingerprint. This is used during resume to ensure the checkpoint
#' is compatible with the current configuration.
#'
#' @param checkpoint A checkpoint list returned by [read_checkpoint()].
#' @param config_fingerprint Character string containing the expected
#'   configuration fingerprint.
#'
#' @return Logical. TRUE if fingerprints match exactly, FALSE otherwise.
#'
#' @details
#' The config fingerprint is a hash of the normalized configuration that
#' excludes runtime-only fields (result_path, checkpoint_every, progress).
#' A mismatch indicates that the checkpoint was created with a different
#' configuration and should not be used for resuming.
#'
#' @seealso [read_checkpoint()]
#'
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' checkpoint <- read_checkpoint("/path/to/results")
#' if (validate_checkpoint_fingerprint(checkpoint, expected_fingerprint)) {
#'   # Safe to resume from this checkpoint
#' } else {
#'   # Configuration has changed, cannot resume
#' }
#' }
validate_checkpoint_fingerprint <- function(checkpoint, config_fingerprint) {
  if (is.null(checkpoint) || is.null(checkpoint$meta)) {
    return(FALSE)
  }
  identical(checkpoint$meta$config_fingerprint, config_fingerprint)
}

# =============================================================================
# Checkpoint Listing
# =============================================================================

#' List available checkpoints
#'
#' Returns a vector of all checkpoint IDs available in the result path.
#'
#' @param result_path Character string giving the base path for results.
#'
#' @return Integer vector of checkpoint IDs, sorted in ascending order.
#'   Returns an empty integer vector if no checkpoints exist or the
#'   checkpoints directory doesn't exist.
#'
#' @details
#' This function scans the checkpoints subdirectory and extracts IDs
#' from directory names matching the pattern `cp_XXXXXX`.
#'
#' @seealso [read_checkpoint()], [get_latest_valid_checkpoint()]
#'
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' checkpoint_ids <- list_checkpoints("/path/to/results")
#' # Returns: c(1L, 2L, 3L, 5L, 10L)
#' }
list_checkpoints <- function(result_path) {
  if (is.null(result_path)) {
    return(integer(0))
  }

  checkpoint_dir <- file.path(result_path, "checkpoints")

  if (!dir.exists(checkpoint_dir)) {
    return(integer(0))
  }

  dirs <- list.dirs(checkpoint_dir, recursive = FALSE, full.names = FALSE)

  # Filter to only valid checkpoint directories (cp_XXXXXX format)
  # This excludes .tmp, etc.
  valid_dirs <- dirs[grepl("^cp_[0-9]{6}$", dirs)]

  if (length(valid_dirs) == 0) {
    return(integer(0))
  }

  # Extract numeric IDs and sort
  checkpoint_ids <- as.integer(gsub("^cp_", "", valid_dirs))

  # Remove NAs from coercion
  checkpoint_ids <- checkpoint_ids[!is.na(checkpoint_ids)]

  sort(checkpoint_ids)
}

# =============================================================================
# Helper: Find Latest Valid Checkpoint
# =============================================================================

#' Get latest valid checkpoint
#'
#' Scans backwards from the latest checkpoint to find the most recent
#' checkpoint with valid checksums. Used for corruption recovery.
#'
#' @param result_path Character string giving the base path for results.
#' @param config_fingerprint Optional character string. If provided, also
#'   validates that the checkpoint's fingerprint matches.
#'
#' @return A checkpoint list (from [read_checkpoint()]), or NULL if no
#'   valid checkpoint is found.
#'
#' @details
#' This function implements the corruption recovery algorithm:
#' \enumerate{
#'   \item Get list of all checkpoint IDs
#'   \item Start from the highest ID
#'   \item Try to read each checkpoint (validates checksums)
#'   \item If fingerprint is provided, also validate fingerprint
#'   \item Return first valid checkpoint found
#'   \item Return NULL if no valid checkpoint found
#' }
#'
#' @seealso [read_checkpoint()], [list_checkpoints()]
#'
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' # Find latest valid checkpoint (any fingerprint)
#' checkpoint <- get_latest_valid_checkpoint("/path/to/results")
#'
#' # Find latest valid checkpoint with matching fingerprint
#' checkpoint <- get_latest_valid_checkpoint(
#'   "/path/to/results",
#'   config_fingerprint = expected_fingerprint
#' )
#' }
get_latest_valid_checkpoint <- function(
  result_path,
  config_fingerprint = NULL
) {
  if (is.null(result_path)) {
    return(NULL)
  }

  checkpoint_ids <- list_checkpoints(result_path)

  if (length(checkpoint_ids) == 0) {
    return(NULL)
  }

  # Try checkpoints from newest to oldest
  for (id in rev(checkpoint_ids)) {
    checkpoint <- read_checkpoint(result_path, checkpoint_id = id)

    if (!is.null(checkpoint)) {
      # If fingerprint validation requested, check it
      if (!is.null(config_fingerprint)) {
        if (!validate_checkpoint_fingerprint(checkpoint, config_fingerprint)) {
          next # Fingerprint mismatch, try older checkpoint
        }
      }
      return(checkpoint)
    }
  }

  NULL
}

# =============================================================================
# Helper: Convert Task Results to Data Frame
# =============================================================================

#' Convert task results to data frame
#'
#' Converts a list of `bayesim_task_result` objects to a data frame
#' suitable for checkpoint storage.
#'
#' @param task_results A list of `bayesim_task_result` objects.
#'
#' @return A data frame with one row per task, containing:
#'   \itemize{
#'     \item `task_id` - Task identifier
#'     \item `status` - Task status ("success", "failed", "skipped")
#'     \item Columns from `metrics` (if present)
#'     \item Columns from `diagnostics` (if present)
#'   }
#'
#' @details
#' This function flattens the nested structure of task results into a
#' tabular format. Metrics and diagnostics are extracted and added as
#' columns. NULL task results are skipped.
#'
#' @keywords internal
#' @export
results_to_dataframe <- function(task_results) {
  if (is.null(task_results) || length(task_results) == 0) {
    return(data.frame(task_id = character(), status = character()))
  }

  if (!is.list(task_results)) {
    stop(bayesim_contract_error("task_results must be a list"))
  }

  rows <- lapply(task_results, function(tr) {
    if (is.null(tr)) {
      return(NULL)
    }

    # Start with basic task info
    row <- list(
      task_id = tr$task_id,
      status = tr$status
    )

    # Add metrics if present (flattening named numeric vectors)
    if (!is.null(tr$metrics) && length(tr$metrics) > 0) {
      row <- c(row, flatten_with_prefix(tr$metrics, ""))
    }

    # Add diagnostics if present (flattening named numeric vectors)
    if (!is.null(tr$diagnostics) && length(tr$diagnostics) > 0) {
      row <- c(row, flatten_with_prefix(tr$diagnostics, ""))
    }

    # Add error info if present
    if (!is.null(tr$error)) {
      row$error_class <- tr$error$error_class
      row$error_message <- tr$error$error_message
    }

    # Add timing if present
    if (!is.null(tr$timing) && !is.null(tr$timing$total)) {
      row$timing_total <- tr$timing$total
    }

    as.data.frame(row, stringsAsFactors = FALSE)
  })

  # Remove NULL entries
  non_null <- vapply(rows, Negate(is.null), logical(1))
  rows <- rows[non_null]

  if (length(rows) == 0) {
    return(data.frame(task_id = character(), status = character()))
  }

  # Use dplyr::bind_rows() which handles different columns gracefully
  # by filling missing columns with NA
  dplyr::bind_rows(rows)
}

# =============================================================================
# Helper: Read Run Manifest
# =============================================================================

#' Read run manifest
#'
#' Reads the run manifest from the result path.
#'
#' @param result_path Character string giving the base path for results.
#'
#' @return A list containing the run manifest, or NULL if not found.
#'
#' @details
#' The run manifest contains:
#' \itemize{
#'   \item `run_schema_version` - Schema version for the run structure
#'   \item `result_schema_version` - Schema version for result files
#'   \item `config_fingerprint` - Hash of the configuration
#'   \item `created` - Timestamp when the run was created
#' }
#'
#' @keywords internal
read_run_manifest <- function(result_path) {
  if (is.null(result_path)) {
    return(NULL)
  }

  manifest_path <- file.path(result_path, "run_manifest.json")

  if (!file.exists(manifest_path)) {
    return(NULL)
  }

  jsonlite::read_json(manifest_path)
}

checkpoint_data_path <- function(directory, stem, checkpoint_format = "rds") {
  if (!identical(checkpoint_format, "rds")) {
    stop(bayesim_checkpoint_error(paste0("Unsupported checkpoint format '", checkpoint_format, "'")))
  }

  file.path(directory, paste0(stem, ".rds"))
}

write_checkpoint_object <- function(x, path, checkpoint_format = "rds") {
  write_rds_atomic(x, path)
}

read_checkpoint_object <- function(path, checkpoint_format = "rds") {
  tryCatch(
    readRDS(path),
    error = function(e) {
      stop(bayesim_checkpoint_error(
        paste0("Failed to read checkpoint file '", path, "': ", conditionMessage(e))
      ))
    }
  )
}

assert_supported_checkpoint_format <- function(checkpoint_format) {
  if (!identical(checkpoint_format, "rds")) {
    stop(bayesim_checkpoint_error(paste0("Unsupported checkpoint format '", checkpoint_format, "'")))
  }

  invisible(checkpoint_format)
}

# =============================================================================
# Helper: Validate Schema Compatibility
# =============================================================================

#' Validate schema compatibility
#'
#' Checks whether the schema versions in a run manifest are compatible
#' with the current package version.
#'
#' @param manifest A run manifest list from [read_run_manifest()].
#'
#' @return Logical. TRUE if schemas are compatible, FALSE otherwise.
#'
#' @details
#' Currently, this function checks for exact version matches. In the future,
#' it may support backward compatibility for certain version ranges.
#'
#' @keywords internal
validate_schema_compatibility <- function(manifest) {
  if (is.null(manifest)) {
    return(FALSE)
  }

  run_version <- manifest$run_schema_version
  result_version <- manifest$result_schema_version

  # Currently require exact version match
  # In the future, this could support backward compatibility
  identical(run_version, RUN_SCHEMA_VERSION) &&
    identical(result_version, RESULT_SCHEMA_VERSION)
}

# =============================================================================
# Helper: Clean Old Checkpoints
# =============================================================================

#' Clean old checkpoints
#'
#' Removes old checkpoints, keeping only the N most recent ones.
#'
#' @param result_path Character string giving the base path for results.
#' @param keep_n Integer. Number of most recent checkpoints to keep.
#'   Default is 5.
#'
#' @return Invisible number of checkpoints removed.
#'
#' @details
#' This function can be used to manage disk space by removing old checkpoints
#' that are no longer needed. It always keeps the latest checkpoint and
#' any checkpoints that are more recent than the keep_n threshold.
#'
#' @keywords internal
#' @export
clean_old_checkpoints <- function(result_path, keep_n = 5) {
  if (is.null(result_path)) {
    return(invisible(0))
  }

  checkpoint_ids <- list_checkpoints(result_path)

  if (length(checkpoint_ids) <= keep_n) {
    return(invisible(0))
  }

  # Identify checkpoints to remove (all but the keep_n most recent)
  ids_to_remove <- sort(checkpoint_ids)[seq_len(
    length(checkpoint_ids) - keep_n
  )]

  removed <- 0
  for (id in ids_to_remove) {
    checkpoint_name <- sprintf("cp_%06d", id)
    checkpoint_dir <- file.path(result_path, "checkpoints", checkpoint_name)

    if (dir.exists(checkpoint_dir)) {
      unlink(checkpoint_dir, recursive = TRUE)
      removed <- removed + 1
    }
  }

  invisible(removed)
}
