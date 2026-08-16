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

# Schema versions are defined in bayesim-package.R
# These are imported via NAMESPACE

# =============================================================================
# Checkpoint Directory Initialization
# =============================================================================

#' Initialize checkpoint directory
#'
#' Creates the checkpoint directory structure for a new simulation run.
#' This includes the base result path, checkpoints and outcomes subdirectories,
#' run manifest, and initial latest pointer. The ledger subdirectory is created
#' lazily at the first delta write.
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
#' +-- checkpoints/         # atomic checkpoint commit directories
#' +-- outcomes/            # immutable, redundantly mirrored outcome shards
#' +-- ledger/              # base ledger plus status deltas (lazy)
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
  checkpoint_format = "rds",
  retention_spec = NULL,
  run_policy_spec = NULL
) {
  if (is.null(result_path)) {
    return(NULL)
  }

  assert_supported_checkpoint_format(checkpoint_format)

  if (file.exists(result_path) && !dir.exists(result_path)) {
    cli::cli_abort(
      "Result path exists but is not a directory: {.file {result_path}}",
      class = "bayesim_checkpoint_error"
    )
  }

  manifest_path <- file.path(result_path, "run_manifest.json")
  if (file.exists(manifest_path)) {
    existing <- tryCatch(
      jsonlite::read_json(manifest_path),
      error = function(e) NULL
    )
    if (is.null(existing)) {
      cli::cli_abort(
        "Result path contains an unreadable run manifest; refusing to overwrite it",
        class = "bayesim_checkpoint_error"
      )
    }
    if (!identical(existing$config_fingerprint, config_fingerprint)) {
      cli::cli_abort(
        c(
          "Result path belongs to a different simulation study",
          x = "Existing and requested study fingerprints differ.",
          i = "Choose a new result_path or explicitly resume the existing study."
        ),
        class = "bayesim_checkpoint_error"
      )
    }
    existing_format <- existing$checkpoint_format %||% "rds"
    if (!identical(existing_format, checkpoint_format)) {
      cli::cli_abort(
        "Checkpoint format mismatch for existing result path",
        class = "bayesim_checkpoint_error"
      )
    }
    dir.create(
      file.path(result_path, "checkpoints"),
      recursive = TRUE,
      showWarnings = FALSE
    )
    dir.create(
      file.path(result_path, "outcomes"),
      recursive = TRUE,
      showWarnings = FALSE
    )
    if (!file.exists(file.path(result_path, "latest.json"))) {
      write_json_atomic(
        list(checkpoint_id = NULL),
        file.path(result_path, "latest.json")
      )
    }
    return(invisible(result_path))
  }

  if (dir.exists(result_path) && !is_empty_result_path(result_path)) {
    cli::cli_abort(
      c(
        "Result path is not empty and does not contain a bayesim run manifest",
        i = "Refusing to mix a new study with existing files."
      ),
      class = "bayesim_checkpoint_error"
    )
  }

  # Create base directory
  dir.create(result_path, recursive = TRUE, showWarnings = FALSE)

  # Create checkpoints subdirectory
  dir.create(
    file.path(result_path, "checkpoints"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  dir.create(
    file.path(result_path, "outcomes"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Write run manifest
  manifest <- list(
    run_schema_version = RUN_SCHEMA_VERSION,
    result_schema_version = RESULT_SCHEMA_VERSION,
    config_fingerprint = config_fingerprint,
    config_spec = config_spec,
    retention_spec = retention_spec,
    run_policy = run_policy_spec,
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
#' @param keep_checkpoints Positive integer; number of checkpoint commit
#'   directories to retain. Older commit directories are pruned only after the
#'   new checkpoint commit and `latest.json` have been written successfully.
#'   Pruning never removes immutable outcome shards or ledger history, so
#'   durable storage grows roughly linearly with completed tasks. `Inf`
#'   disables pruning for direct/internal callers.
#' @param prior_results_df Optional cached data frame of results from before the
#'   current execution. When supplied, it replaces the legacy read of the
#'   previous checkpoint.
#' @param prior_task_results Optional canonical task outcomes from before the
#'   current execution, used when migrating legacy checkpoints.
#' @param adaptive_next_check Optional persisted adaptive-check threshold.
#' @param adaptive_state Optional persistable adaptive precision snapshot.
#' @param run_policy_spec Optional serialized effective run policy.
#' @param prior_checkpoint Optional validated in-memory checkpoint state. The
#'   filesystem store supplies this to keep live append work proportional to
#'   newly completed tasks.
#' @param return_state Logical; return the validated checkpoint state instead
#'   of only its numeric ID.
#' @param delta_store Logical; write immutable outcome and ledger deltas instead
#'   of legacy full snapshots.
#'
#' @return Invisible checkpoint ID (integer), or NULL if result_path is NULL.
#'
#' @details
#' The checkpoint directory structure is:
#' \preformatted{
#' checkpoints/
#' +-- cp_000001/
#'     +-- meta.json         # checkpoint metadata
#'     +-- ledger.rds        # ledger view (delta in delta-store mode)
#'     +-- results.rds       # results view (delta in delta-store mode)
#'     +-- checksums.json    # file integrity checksums
#' }
#'
#' With `delta_store = TRUE`, each checkpoint commit directory is an atomic
#' commit record: its `meta.json` plus the delta views above. Newly completed
#' outcomes are appended once as immutable, redundantly mirrored shards under
#' `outcomes/`, and task statuses under `ledger/` are a base ledger plus
#' status deltas. Readers reconstruct the accumulated run state from the
#' shards a commit references; they never rewrite history. With
#' `delta_store = FALSE` (legacy snapshot mode for direct/internal callers),
#' the checkpoint directory itself carries the snapshot views.
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
  checkpoint_format = "rds",
  keep_checkpoints = Inf,
  prior_results_df = NULL,
  prior_task_results = NULL,
  adaptive_next_check = NULL,
  adaptive_state = NULL,
  run_policy_spec = NULL,
  prior_checkpoint = NULL,
  return_state = FALSE,
  delta_store = FALSE
) {
  if (is.null(result_path)) {
    return(invisible(NULL))
  }

  assert_supported_checkpoint_format(checkpoint_format)

  has_supplied_previous <- !is.null(prior_checkpoint)
  previous <- prior_checkpoint
  if (is.null(previous) && can_resume(result_path)) {
    previous <- get_latest_valid_checkpoint(
      result_path,
      config_fingerprint = config_fingerprint,
      load_outcomes = FALSE
    )
  }

  current_outcomes <- Filter(
    function(x) !is.null(x) && is_bayesim_task_result(x),
    task_results
  )
  previous_grid <- previous$task_grid %||%
    data.frame(
      task_id = character(),
      status = character()
    )
  previous_status <- stats::setNames(
    previous_grid$status,
    previous_grid$task_id
  )
  previous_reason <- if ("stop_reason" %in% names(previous_grid)) {
    stats::setNames(previous_grid$stop_reason, previous_grid$task_id)
  } else {
    stats::setNames(
      rep(NA_character_, nrow(previous_grid)),
      previous_grid$task_id
    )
  }
  changed <- vapply(
    current_outcomes,
    function(outcome) {
      if (!has_supplied_previous) {
        return(TRUE)
      }
      old_status <- unname(previous_status[outcome$task_id])
      old_reason <- unname(previous_reason[outcome$task_id])
      new_reason <- outcome$stop_reason %||% NA_character_
      length(old_status) == 0L ||
        is.na(old_status) ||
        !identical(old_status, outcome$status) ||
        !identical(old_reason, new_reason)
    },
    logical(1)
  )
  new_outcomes <- current_outcomes[changed]

  # The final assembly often asks the store to write the exact ledger already
  # checkpointed by the last execution batch. Avoid a redundant snapshot: it
  # would add no durability and could prune the predecessor needed for corrupt
  # latest-shard fallback.
  same_ledger <- !is.null(previous) &&
    identical(
      as.data.frame(previous_grid),
      as.data.frame(task_grid)
    )
  same_policy <- identical(
    previous$meta$run_policy %||% NULL,
    run_policy_spec %||% NULL
  )
  if (
    has_supplied_previous &&
      same_policy &&
      same_ledger &&
      length(new_outcomes) == 0L
  ) {
    return(invisible(
      if (isTRUE(return_state)) {
        previous
      } else {
        previous$checkpoint_id
      }
    ))
  }

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
  new_shard_path <- NULL
  new_ledger_path <- NULL

  # Use tryCatch to ensure cleanup on error
  tryCatch(
    {
      outcome_refs <- normalize_outcome_shard_refs(
        previous$meta$outcome_shards %||% NULL
      )

      # One-time migration for checkpoints written before append-only outcome
      # shards existed. Prefer canonical outcomes when available; otherwise
      # reconstruct the legacy flat rows while keeping truth/error fields out
      # of metrics.
      if (
        !isTRUE(delta_store) &&
          !is.null(previous) &&
          length(outcome_refs) == 0L &&
          length(current_outcomes) > 0L
      ) {
        legacy_outcomes <- previous$task_outcomes %||% prior_task_results
        if (is.null(legacy_outcomes) || length(legacy_outcomes) == 0L) {
          legacy_outcomes <- task_outcomes_from_dataframe(
            previous$results_df %||% prior_results_df
          )
        }
        legacy_outcomes <- Filter(
          function(x) !is.null(x) && is_bayesim_task_result(x),
          legacy_outcomes %||% list()
        )
        if (length(legacy_outcomes) > 0L) {
          new_ids <- vapply(new_outcomes, function(x) x$task_id, character(1))
          legacy_outcomes <- legacy_outcomes[
            !vapply(
              legacy_outcomes,
              function(x) x$task_id %in% new_ids,
              logical(1)
            )
          ]
          new_outcomes <- c(legacy_outcomes, new_outcomes)
        }
      }

      if (length(new_outcomes) > 0L) {
        shard_name <- sprintf("shard_%06d.rds", checkpoint_id)
        shard_path <- file.path(result_path, "outcomes", shard_name)
        new_shard_path <- shard_path
        new_ref <- if (isTRUE(delta_store)) {
          write_redundant_shard(new_outcomes, shard_path)
        } else {
          write_rds_atomic(new_outcomes, shard_path)
          list(
            file = shard_name,
            checksum = compute_file_checksum(shard_path),
            n_outcomes = length(new_outcomes)
          )
        }
        new_ref$n_outcomes <- length(new_outcomes)
        outcome_refs <- c(outcome_refs, list(new_ref))
      }

      ledger_to_write <- task_grid
      ledger_ref <- NULL
      if (isTRUE(delta_store)) {
        ledger_dir <- file.path(result_path, "ledger")
        dir.create(ledger_dir, recursive = TRUE, showWarnings = FALSE)
        if (is.null(previous)) {
          ledger_name <- "base_000001.rds"
          ledger_to_write <- task_grid
        } else {
          ledger_name <- sprintf("delta_%06d.rds", checkpoint_id)
          old_status <- stats::setNames(
            previous_grid$status,
            previous_grid$task_id
          )
          old_reason <- if ("stop_reason" %in% names(previous_grid)) {
            stats::setNames(previous_grid$stop_reason, previous_grid$task_id)
          } else {
            stats::setNames(
              rep(NA_character_, nrow(previous_grid)),
              previous_grid$task_id
            )
          }
          status_changed <- vapply(
            seq_len(nrow(task_grid)),
            function(i) {
              id <- task_grid$task_id[[i]]
              reason <- if ("stop_reason" %in% names(task_grid)) {
                task_grid$stop_reason[[i]]
              } else {
                NA_character_
              }
              !identical(unname(old_status[id]), task_grid$status[[i]]) ||
                !identical(unname(old_reason[id]), reason)
            },
            logical(1)
          )
          ledger_to_write <- task_grid[
            status_changed,
            intersect(
              c("task_id", "status", "stop_reason"),
              names(task_grid)
            ),
            drop = FALSE
          ]
        }
        ledger_ref <- write_redundant_shard(
          ledger_to_write,
          file.path(ledger_dir, ledger_name)
        )
        new_ledger_path <- file.path(ledger_dir, ledger_name)
        ledger_ref$n_rows <- nrow(ledger_to_write)
      }

      # Write meta.json
      meta <- list(
        checkpoint_id = checkpoint_id,
        created = as.character(Sys.time()),
        config_fingerprint = config_fingerprint,
        n_tasks = nrow(task_grid),
        n_success = sum(task_grid$status == "success", na.rm = TRUE),
        n_failed = sum(task_grid$status == "failed", na.rm = TRUE),
        n_pending = sum(
          is_resumable_task_status(task_grid$status),
          na.rm = TRUE
        ),
        n_policy_stopped = if ("stop_reason" %in% names(task_grid)) {
          sum(
            task_grid$status == "skipped" &
              task_grid$stop_reason %in% c("max_errors", "adaptive_stop"),
            na.rm = TRUE
          )
        } else {
          0L
        },
        adaptive_next_check = adaptive_next_check,
        adaptive_state = adaptive_state,
        run_policy = run_policy_spec,
        storage_mode = if (isTRUE(delta_store)) "delta-v1" else "snapshot-v1",
        ledger_shard = ledger_ref,
        outcome_shard = if (isTRUE(delta_store) && length(new_outcomes) > 0L) {
          utils::tail(outcome_refs, 1L)[[1L]]
        } else {
          NULL
        },
        outcome_shards = if (isTRUE(delta_store)) NULL else outcome_refs
      )
      write_json_atomic(meta, file.path(tmp_dir, "meta.json"))

      # Only the "rds" format is supported (asserted above), so the data files
      # are addressed and (de)serialized directly.
      ledger_path <- file.path(tmp_dir, "ledger.rds")
      results_path <- file.path(tmp_dir, "results.rds")

      # Write ledger (task grid)
      write_rds_atomic(ledger_to_write, ledger_path)

      # The flat checkpoint view is a delta matching this checkpoint's new
      # canonical shard. read_checkpoint() derives the accumulated flat view
      # from all referenced shards, so write cost stays proportional to new
      # outcomes rather than total study size.
      results_df <- if (length(new_outcomes) > 0L) {
        results_to_dataframe(new_outcomes)
      } else {
        # Compatibility for direct/internal callers that still supply the
        # historical unclassed task-result lists. The execution engine always
        # supplies canonical outcomes and therefore takes the delta path.
        results_to_dataframe(task_results)
      }
      write_rds_atomic(results_df, results_path)

      # Read-back validation
      test_grid <- tryCatch(readRDS(ledger_path), error = function(e) NULL)
      if (is.null(test_grid)) {
        stop(bayesim_checkpoint_error(
          "Checkpoint read-back validation failed for ledger"
        ))
      }

      test_results <- tryCatch(readRDS(results_path), error = function(e) NULL)
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

      prune_checkpoints(result_path, keep = keep_checkpoints)

      checkpoint_state <- list(
        checkpoint_id = checkpoint_id,
        meta = meta,
        task_grid = task_grid,
        results_df = results_df,
        task_outcomes = NULL,
        adaptive_next_check = adaptive_next_check,
        adaptive_state = adaptive_state,
        run_policy = run_policy_spec
      )
      invisible(if (isTRUE(return_state)) checkpoint_state else checkpoint_id)
    },
    error = function(e) {
      # Clean up temporary directory on any error
      if (dir.exists(tmp_dir)) {
        unlink(tmp_dir, recursive = TRUE)
      }
      # A failed commit must leave no half-published files outside the commit
      # directory either: remove the outcome shard and the ledger shard this
      # attempt wrote, together with their redundant mirrors and checksum
      # descriptors (.mirror.rds / .json / .mirror.json siblings). Only do so
      # when the commit directory itself never materialized — otherwise the
      # files belong to an existing checkpoint commit.
      if (!dir.exists(checkpoint_dir)) {
        shard_stems <- character()
        if (!is.null(new_shard_path)) {
          shard_stems <- c(shard_stems, sub("\\.rds$", "", new_shard_path))
        }
        if (!is.null(new_ledger_path)) {
          shard_stems <- c(shard_stems, sub("\\.rds$", "", new_ledger_path))
        }
        for (stem in unique(shard_stems)) {
          unlink(paste0(
            stem,
            c(".rds", ".mirror.rds", ".json", ".mirror.json")
          ))
        }
      }
      stop(e)
    }
  )
}

#' Prune old checkpoint commit directories
#'
#' @param result_path Checkpoint result directory.
#' @param keep Number of newest checkpoint commits to keep; `Inf` keeps all.
#'   Removes commit directories only; immutable outcome shards and ledger
#'   history are never pruned.
#' @return Invisible vector of removed checkpoint IDs.
#' @keywords internal
prune_checkpoints <- function(result_path, keep = 2L) {
  if (is.infinite(keep)) {
    return(invisible(integer()))
  }
  keep <- as.integer(keep)
  if (length(keep) != 1L || is.na(keep) || keep < 1L) {
    stop(bayesim_config_error("keep_checkpoints must be a positive integer"))
  }

  checkpoint_ids <- list_checkpoints(result_path)
  if (length(checkpoint_ids) <= keep) {
    return(invisible(integer()))
  }

  remove_ids <- utils::head(checkpoint_ids, -keep)
  remove_paths <- file.path(
    result_path,
    "checkpoints",
    sprintf("cp_%06d", remove_ids)
  )
  unlink(remove_paths, recursive = TRUE)
  invisible(remove_ids)
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
#' @param load_outcomes Logical; when FALSE, validate shard integrity but omit
#'   materializing accumulated outcomes. Used by the filesystem RunStore while
#'   appending the next ledger checkpoint.
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
#' @examples
#' \dontrun{
#' # Read latest checkpoint
#' checkpoint <- read_checkpoint("/path/to/results")
#'
#' # Read specific checkpoint
#' checkpoint <- read_checkpoint("/path/to/results", checkpoint_id = 5)
#' }
read_checkpoint <- function(
  result_path,
  checkpoint_id = NULL,
  load_outcomes = TRUE
) {
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

  # Read checkpoint files (the supported format is always "rds", asserted via
  # the manifest above, so the files are addressed directly)
  meta <- jsonlite::read_json(file.path(checkpoint_dir, "meta.json"))
  task_grid <- readRDS(file.path(checkpoint_dir, "ledger.rds"))
  checkpoint_results_df <- readRDS(file.path(checkpoint_dir, "results.rds"))

  delta_store <- identical(meta$storage_mode %||% NULL, "delta-v1")
  if (delta_store) {
    task_grid <- read_delta_ledger(result_path, checkpoint_id)
    if (is.null(task_grid)) {
      cli::cli_warn(
        "Checkpoint {checkpoint_id} references a missing or corrupt ledger shard"
      )
      return(NULL)
    }
  }

  outcome_refs <- if (delta_store) {
    list_delta_shard_refs(
      file.path(result_path, "outcomes"),
      "shard_",
      checkpoint_id
    )
  } else {
    normalize_outcome_shard_refs(meta$outcome_shards %||% NULL)
  }
  if (length(outcome_refs) > 0L) {
    shard_read <- read_outcome_shards(
      result_path,
      outcome_refs,
      load_outcomes = load_outcomes
    )
    if (!isTRUE(shard_read$valid)) {
      cli::cli_warn(
        "Checkpoint {checkpoint_id} references a missing or corrupt outcome shard"
      )
      return(NULL)
    }
    task_outcomes <- shard_read$outcomes
    results_df <- if (isTRUE(load_outcomes)) {
      results_to_dataframe(task_outcomes)
    } else {
      checkpoint_results_df
    }
  } else {
    # Legacy schema: full flat results and optionally full structured outcomes
    # lived inside each checkpoint directory.
    outcomes_path <- file.path(checkpoint_dir, "outcomes.rds")
    task_outcomes <- if (isTRUE(load_outcomes) && file.exists(outcomes_path)) {
      readRDS(outcomes_path)
    } else {
      NULL
    }
    results_df <- checkpoint_results_df
  }

  list(
    checkpoint_id = checkpoint_id,
    meta = meta,
    task_grid = task_grid,
    results_df = results_df,
    task_outcomes = task_outcomes,
    adaptive_next_check = meta$adaptive_next_check %||% NULL,
    adaptive_state = meta$adaptive_state %||% NULL
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
  config_fingerprint = NULL,
  load_outcomes = TRUE
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
    checkpoint <- read_checkpoint(
      result_path,
      checkpoint_id = id,
      load_outcomes = load_outcomes
    )

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

#' Normalize outcome-shard records read from JSON metadata.
#' @keywords internal
normalize_outcome_shard_refs <- function(refs) {
  if (is.null(refs) || length(refs) == 0L) {
    return(list())
  }
  if (is.character(refs)) {
    return(lapply(refs, function(file) list(file = file, checksum = NULL)))
  }
  if (is.list(refs) && !is.null(refs$file)) {
    refs <- list(refs)
  }
  lapply(refs, function(ref) {
    if (is.character(ref)) {
      return(list(file = ref, checksum = NULL))
    }
    list(
      file = ref$file %||% NA_character_,
      checksum = ref$checksum %||% NULL,
      mirror = ref$mirror %||% NULL,
      mirror_checksum = ref$mirror_checksum %||% NULL,
      n_outcomes = ref$n_outcomes %||% NULL
    )
  })
}

#' Write one immutable shard with an independently checksummed mirror.
#'
#' The descriptor is committed last, so readers never observe a partially
#' published shard. Mirroring keeps later checkpoints readable after damage to
#' any single historical shard without rewriting accumulated history.
#' @keywords internal
write_redundant_shard <- function(object, path) {
  mirror <- sub("\\.rds$", ".mirror.rds", path)
  descriptor <- sub("\\.rds$", ".json", path)
  descriptor_mirror <- sub("\\.rds$", ".mirror.json", path)
  write_rds_atomic(object, path)
  write_rds_atomic(object, mirror)
  ref <- list(
    file = basename(path),
    checksum = compute_file_checksum(path),
    mirror = basename(mirror),
    mirror_checksum = compute_file_checksum(mirror)
  )
  write_json_atomic(ref, descriptor)
  write_json_atomic(ref, descriptor_mirror)
  ref
}

#' Read a checksummed shard, falling back to its mirror.
#' @keywords internal
read_redundant_shard <- function(directory, ref) {
  candidates <- list(
    list(file = ref$file, checksum = ref$checksum),
    list(file = ref$mirror, checksum = ref$mirror_checksum)
  )
  for (i in seq_along(candidates)) {
    candidate <- candidates[[i]]
    if (is.null(candidate$file) || is.null(candidate$checksum)) {
      next
    }
    path <- file.path(directory, candidate$file)
    valid <- file.exists(path) &&
      identical(
        compute_file_checksum(path),
        candidate$checksum
      )
    if (isTRUE(valid)) {
      value <- tryCatch(readRDS(path), error = function(e) NULL)
      if (!is.null(value)) {
        if (i == 2L) {
          cli::cli_warn("Recovered a corrupt shard from its mirror: {ref$file}")
        }
        return(list(valid = TRUE, value = value))
      }
    }
  }
  list(valid = FALSE, value = NULL)
}

#' List committed delta-store shard descriptors through a checkpoint.
#' @keywords internal
list_delta_shard_refs <- function(directory, prefix, checkpoint_id) {
  if (!dir.exists(directory)) {
    return(list())
  }
  files <- list.files(
    directory,
    pattern = paste0(
      "^",
      prefix,
      "[0-9]{6}(\\.rds|\\.mirror\\.rds|\\.json|\\.mirror\\.json)$"
    ),
    full.names = FALSE
  )
  ids <- as.integer(sub(
    paste0("^", prefix, "([0-9]{6}).*$"),
    "\\1",
    files
  ))
  ids <- sort(unique(ids[!is.na(ids) & ids <= checkpoint_id]))
  refs <- lapply(ids, function(id) {
    stem <- sprintf("%s%06d", prefix, id)
    descriptors <- file.path(
      directory,
      c(paste0(stem, ".json"), paste0(stem, ".mirror.json"))
    )
    for (i in seq_along(descriptors)) {
      path <- descriptors[[i]]
      ref <- tryCatch(jsonlite::read_json(path), error = function(e) NULL)
      if (!is.null(ref)) {
        if (i == 2L) {
          cli::cli_warn(
            "Recovered a corrupt shard descriptor from its mirror: {stem}"
          )
        }
        return(ref)
      }
    }
    # Preserve the missing/corrupt shard in the returned sequence so callers
    # fail closed instead of silently assembling an incomplete run.
    list(file = NA_character_, checksum = NULL)
  })
  refs
}

#' Reconstruct the task ledger from its immutable base and status deltas.
#' @keywords internal
read_delta_ledger <- function(result_path, checkpoint_id) {
  ledger_dir <- file.path(result_path, "ledger")
  base_refs <- list_delta_shard_refs(ledger_dir, "base_", checkpoint_id)
  if (length(base_refs) != 1L) {
    return(NULL)
  }
  base <- read_redundant_shard(ledger_dir, base_refs[[1L]])
  if (!isTRUE(base$valid)) {
    return(NULL)
  }
  ledger <- base$value
  deltas <- list_delta_shard_refs(ledger_dir, "delta_", checkpoint_id)
  for (ref in deltas) {
    delta <- read_redundant_shard(ledger_dir, ref)
    if (!isTRUE(delta$valid)) {
      return(NULL)
    }
    if (nrow(delta$value) == 0L) {
      next
    }
    hit <- match(delta$value$task_id, ledger$task_id)
    if (anyNA(hit)) {
      return(NULL)
    }
    for (field in intersect(c("status", "stop_reason"), names(delta$value))) {
      ledger[[field]][hit] <- delta$value[[field]]
    }
  }
  ledger
}

#' Validate and optionally materialize accumulated canonical outcome shards.
#' @keywords internal
read_outcome_shards <- function(
  result_path,
  refs,
  load_outcomes = TRUE
) {
  refs <- normalize_outcome_shard_refs(refs)
  outcomes <- list()
  outcome_ids <- character()
  for (ref in refs) {
    file <- ref$file
    if (
      !is.character(file) ||
        length(file) != 1L ||
        is.na(file) ||
        !identical(basename(file), file)
    ) {
      return(list(valid = FALSE, outcomes = NULL))
    }
    if (!is.null(ref$mirror)) {
      recovered <- read_redundant_shard(
        file.path(result_path, "outcomes"),
        ref
      )
      if (!isTRUE(recovered$valid)) {
        return(list(valid = FALSE, outcomes = NULL))
      }
      shard <- recovered$value
    } else {
      path <- file.path(result_path, "outcomes", file)
      if (!file.exists(path)) {
        return(list(valid = FALSE, outcomes = NULL))
      }
      checksum <- ref$checksum
      if (
        !is.null(checksum) &&
          !identical(compute_file_checksum(path), checksum)
      ) {
        return(list(valid = FALSE, outcomes = NULL))
      }
      shard <- if (isTRUE(load_outcomes)) {
        tryCatch(readRDS(path), error = function(e) NULL)
      } else {
        list()
      }
    }
    if (!isTRUE(load_outcomes)) {
      next
    }
    if (
      is.null(shard) ||
        !all(vapply(shard, is_bayesim_task_result, logical(1)))
    ) {
      return(list(valid = FALSE, outcomes = NULL))
    }
    # Shards are ordered oldest to newest. Later outcomes replace policy-stop
    # placeholders for the same task identity while preserving deterministic
    # task order in the accumulated view.
    for (outcome in shard) {
      hit <- match(outcome$task_id, outcome_ids)
      if (is.na(hit)) {
        outcomes <- c(outcomes, list(outcome))
        outcome_ids <- c(outcome_ids, outcome$task_id)
      } else {
        outcomes[[hit]] <- outcome
      }
    }
  }
  list(valid = TRUE, outcomes = if (isTRUE(load_outcomes)) outcomes else NULL)
}

#' Migrate the legacy flat checkpoint view to canonical task outcomes.
#'
#' The flat view produced by [results_to_dataframe()] lays columns out in a
#' fixed order: `task_id`, `status`, `stop_reason`, metric columns,
#' `truth__*` columns, diagnostic columns, then the optional
#' `error_class`/`error_message`/`timing_total` trailer. Metric and diagnostic
#' columns share a namespace, so the two groups are separated positionally
#' around the `truth__` block. When a row carries no truth block the groups are
#' contiguous and indistinguishable; every value column is then restored as a
#' metric so the flat representation still round-trips unchanged.
#'
#' @keywords internal
task_outcomes_from_dataframe <- function(results_df) {
  if (is.null(results_df) || nrow(results_df) == 0L) {
    return(list())
  }
  lapply(seq_len(nrow(results_df)), function(i) {
    row <- results_df[i, , drop = FALSE]
    status <- as.character(row$status[[1]])
    truth_cols <- grep("^truth__", names(row), value = TRUE)
    truth <- if (length(truth_cols) > 0L) {
      value <- unlist(row[truth_cols], use.names = FALSE)
      names(value) <- sub("^truth__", "", truth_cols)
      value
    } else {
      NULL
    }
    excluded <- c(
      "task_id",
      "status",
      "stop_reason",
      "error_class",
      "error_message",
      "timing_total",
      truth_cols
    )
    value_cols <- setdiff(names(row), excluded)
    metric_cols <- value_cols
    diagnostic_cols <- character()
    if (length(truth_cols) > 0L) {
      value_positions <- match(value_cols, names(row))
      truth_positions <- match(truth_cols, names(row))
      metric_cols <- value_cols[value_positions < min(truth_positions)]
      diagnostic_cols <- value_cols[value_positions > max(truth_positions)]
    }
    restore <- function(cols) {
      if (length(cols) == 0L) {
        NULL
      } else {
        stats::setNames(lapply(cols, function(name) row[[name]][[1]]), cols)
      }
    }
    metrics <- if (identical(status, "success")) restore(metric_cols) else NULL
    diagnostics <- if (identical(status, "success")) {
      restore(diagnostic_cols)
    } else {
      NULL
    }
    stop_reason <- if (
      "stop_reason" %in% names(row) && !is.na(row$stop_reason[[1]])
    ) {
      as.character(row$stop_reason[[1]])
    } else {
      NULL
    }
    has_error <- "error_class" %in% names(row) && !is.na(row$error_class[[1]])
    error <- if (has_error || identical(status, "failed")) {
      list(
        error_class = if (has_error) {
          as.character(row$error_class[[1]])
        } else {
          "unknown"
        },
        error_message = if (
          "error_message" %in% names(row) && !is.na(row$error_message[[1]])
        ) {
          as.character(row$error_message[[1]])
        } else {
          "Task failed"
        }
      )
    } else {
      NULL
    }
    timing <- if (
      "timing_total" %in% names(row) && is.finite(row$timing_total[[1]])
    ) {
      as.numeric(row$timing_total[[1]])
    } else {
      0
    }
    new_task_result(
      task_id = as.character(row$task_id[[1]]),
      status = status,
      metrics = metrics,
      diagnostics = diagnostics,
      timing = list(total = timing),
      error = error,
      truth = truth,
      stop_reason = stop_reason
    )
  })
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
results_to_dataframe <- function(task_results) {
  if (is.null(task_results) || length(task_results) == 0) {
    return(data.frame(task_id = character(), status = character()))
  }

  rows <- lapply(task_results, function(tr) {
    if (is.null(tr)) {
      return(NULL)
    }

    # Start with basic task info
    row <- list(
      task_id = tr$task_id,
      status = tr$status,
      stop_reason = tr$stop_reason %||% NA_character_
    )

    # Add metrics if present (flattening named numeric vectors)
    if (!is.null(tr$metrics) && length(tr$metrics) > 0) {
      row <- c(row, flatten_with_prefix(tr$metrics, ""))
    }

    # E1: data-generating truth, flattened as truth__<param> columns (always
    # retained — tiny, enables parameter-recovery analysis and plot_recovery()).
    if (!is.null(tr$truth) && length(tr$truth) > 0) {
      tp <- tr$truth
      nm <- if (!is.null(names(tp))) {
        names(tp)
      } else {
        paste0("param", seq_along(tp))
      }
      truth_cols <- as.list(unname(tp))
      names(truth_cols) <- paste0("truth__", nm)
      row <- c(row, truth_cols)
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

  # Combine rows: base-R row-binding that fills missing columns with NA,
  # matching dplyr::bind_rows() semantics without the dplyr dependency.
  all_cols <- unique(unlist(lapply(rows, names)))
  for (i in seq_along(rows)) {
    missing <- setdiff(all_cols, names(rows[[i]]))
    if (length(missing) > 0L) {
      rows[[i]][missing] <- rep(NA, length(missing))
    }
    rows[[i]] <- rows[[i]][all_cols]
  }
  do.call(rbind, rows)
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

assert_supported_checkpoint_format <- function(checkpoint_format) {
  if (!identical(checkpoint_format, "rds")) {
    cli::cli_abort("Unsupported checkpoint format '{checkpoint_format}'")
  }

  invisible(checkpoint_format)
}
