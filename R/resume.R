#' @keywords internal
NULL

#' Check if Resumable Run Exists
#'
#' Checks whether a valid checkpoint exists that can be resumed.
#' A run can be resumed if:
#' - Both run_manifest.json and latest.json exist
#' - latest.json points to a valid checkpoint_id
#' - The referenced checkpoint can be read and validated
#'
#' This is a cheap existence/validity probe: the checkpoint is validated
#' (checksums, ledger, shard integrity) with `load_outcomes = FALSE`, so the
#' full outcome history is never deserialized here. Callers that need the
#' outcomes use [load_for_resume()] instead.
#'
#' @param result_path Character; path to results directory containing checkpoints.
#'
#' @return TRUE if a valid run can be resumed, FALSE otherwise.
#'
#' @keywords internal
can_resume <- function(result_path) {
  if (is.null(result_path)) {
    return(FALSE)
  }

  manifest_path <- file.path(result_path, "run_manifest.json")
  latest_path <- file.path(result_path, "latest.json")

  if (!file.exists(manifest_path) || !file.exists(latest_path)) {
    return(FALSE)
  }

  # Check for valid checkpoint reference
  latest <- tryCatch(
    jsonlite::read_json(latest_path),
    error = function(e) NULL
  )

  if (is.null(latest) || is.null(latest$checkpoint_id)) {
    return(FALSE)
  }

  # Verify checkpoint can be read (without materializing outcomes)
  checkpoint <- tryCatch(
    read_checkpoint(result_path, latest$checkpoint_id, load_outcomes = FALSE),
    error = function(e) NULL
  )

  !is.null(checkpoint)
}

#' Load Run for Resume
#'
#' Loads previous run state for resumption. Validates schema compatibility
#' and configuration fingerprint before allowing resume.
#'
#' @param result_path Character; path to results directory containing checkpoints.
#' @param config SimulationConfig; current configuration object.
#' @param run_store Optional internal [new_run_store()] adapter. When supplied,
#'   it owns checkpoint reads and corruption fallback.
#'
#' @return A list with elements:
#'   - `task_grid`: Task grid with restored status from checkpoint
#'   - `prior_results`: Data frame of results from checkpoint
#'   - `checkpoint_id`: ID of the checkpoint being resumed from
#'
#' @details
#' The function performs the following validation steps:
#' 1. Reads and validates run_manifest.json
#' 2. Checks schema version compatibility
#' 3. Computes and compares configuration fingerprints
#' 4. Finds the most recent valid checkpoint
#' 5. Rebuilds task grid with status from checkpoint
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' config <- simulation_config(...)
#' resume_state <- load_for_resume("/path/to/results", config)
#' task_grid <- resume_state$task_grid
#' prior_results <- resume_state$prior_results
#' }
load_for_resume <- function(result_path, config, run_store = NULL) {
  # Read manifest
  manifest_path <- file.path(result_path, "run_manifest.json")
  manifest <- tryCatch(
    jsonlite::read_json(manifest_path),
    error = function(e) {
      cli::cli_abort(
        "Cannot read run manifest: {manifest_path}",
        class = "bayesim_checkpoint_error"
      )
    }
  )

  # Check run schema version compatibility
  if (!identical(manifest$run_schema_version, RUN_SCHEMA_VERSION)) {
    cli::cli_abort(
      c(
        "Incompatible checkpoint schema",
        "x" = "Checkpoint run_schema_version: {manifest$run_schema_version}",
        "i" = "Expected version: {RUN_SCHEMA_VERSION}"
      ),
      class = "bayesim_checkpoint_error"
    )
  }

  # Check result schema version compatibility
  if (!identical(manifest$result_schema_version, RESULT_SCHEMA_VERSION)) {
    cli::cli_abort(
      c(
        "Incompatible result schema",
        "x" = "Checkpoint result_schema_version: {manifest$result_schema_version}",
        "i" = "Expected version: {RESULT_SCHEMA_VERSION}"
      ),
      class = "bayesim_checkpoint_error"
    )
  }

  manifest_format <- manifest$checkpoint_format %||% "rds"
  if (!identical(manifest_format, config@checkpoint_format)) {
    cli::cli_abort(
      c(
        "Checkpoint format mismatch",
        "x" = "Checkpoint uses '{manifest_format}' but config requests '{config@checkpoint_format}'"
      ),
      class = "bayesim_checkpoint_error"
    )
  }

  # Find valid checkpoint (scan backward from latest if corrupted)
  expected_fingerprint <- compute_config_fingerprint(new_study_spec(config))
  checkpoint <- if (!is.null(run_store)) {
    if (!is_run_store(run_store)) {
      cli::cli_abort("run_store must be a bayesim RunStore adapter")
    }
    run_store$read()
  } else {
    get_latest_valid_checkpoint(
      result_path,
      config_fingerprint = expected_fingerprint
    )
  }
  if (is.null(checkpoint)) {
    cli::cli_abort(
      "No valid compatible checkpoint found in {result_path}",
      class = "bayesim_checkpoint_error"
    )
  }

  validate_resume_retention(
    requested = config@retain,
    persisted = checkpoint$meta$run_policy$retain %||%
      manifest$run_policy$retain %||%
      manifest$retention_spec,
    checkpoint = checkpoint
  )

  # Rebuild task grid with restored status
  fresh_grid <- create_task_grid(config)

  # Merge status from checkpoint ledger
  task_grid <- merge_task_grid_status(fresh_grid, checkpoint$task_grid)

  list(
    task_grid = task_grid,
    prior_results = checkpoint$results_df,
    prior_task_results = checkpoint$task_outcomes,
    adaptive_next_check = checkpoint$adaptive_next_check,
    adaptive_state = checkpoint$adaptive_state,
    checkpoint_id = checkpoint$checkpoint_id
  )
}

#' Validate that resume does not promise artifacts that were never persisted.
#'
#' Retention is runtime policy and therefore does not alter study identity, but
#' widening it after terminal outcomes exist cannot recreate discarded draws,
#' fits, predictions, data, diagnostics, or warnings. Reject that ambiguity at
#' the resume seam before any pending task is executed.
#' @keywords internal
validate_resume_retention <- function(requested, persisted, checkpoint) {
  terminal <- is_terminal_task_status(checkpoint$task_grid$status)
  if (!any(terminal, na.rm = TRUE)) {
    return(invisible(TRUE))
  }

  requested <- resolve_retention_spec(requested)
  if (is.null(persisted) || length(persisted) == 0L) {
    # Legacy manifests did not record retention policy. The historical default
    # guaranteed metrics and diagnostics; heavy artifacts may already have
    # been discarded and therefore cannot safely be promised on resume.
    persisted <- list(
      success = c("metrics", "diagnostics"),
      warning = c("metrics", "diagnostics"),
      error = c("metrics", "diagnostics")
    )
  } else {
    persisted <- resolve_retention_spec(normalize_manifest_retain(persisted))
  }

  widened <- unique(unlist(
    lapply(RETENTION_CONTEXTS, function(context) {
      setdiff(requested[[context]], persisted[[context]])
    }),
    use.names = FALSE
  ))
  widened <- setdiff(widened, "metrics")
  if (length(widened) > 0L) {
    cli::cli_abort(
      c(
        "Cannot widen retention while resuming completed outcomes",
        x = "The existing run did not retain: {paste(widened, collapse = ', ')}.",
        i = "Use the original retention policy or start a new result_path."
      ),
      class = "bayesim_checkpoint_error"
    )
  }
  invisible(TRUE)
}

#' Merge Task Grid Status from Checkpoint
#'
#' Updates a freshly-generated task grid with terminal status values
#' from a checkpoint. Pending tasks remain pending; only tasks that
#' were completed (success/failed/skipped) in the checkpoint are updated.
#'
#' @param fresh_grid Task grid tibble from [create_task_grid()].
#' @param checkpoint_grid Task grid tibble from checkpoint.
#'
#' @return Task grid tibble with status merged from checkpoint.
#'
#' @details
#' The function preserves the deterministic task ordering and RNG seeds
#' from the fresh grid, only updating status for tasks that were terminal
#' in the checkpoint. This ensures resumed runs maintain identical
#' reproducibility guarantees.
#'
#' @keywords internal
merge_task_grid_status <- function(fresh_grid, checkpoint_grid) {
  # Start with fresh grid (has all tasks with pending status and rng_seed)
  # Update status for tasks that were terminal in checkpoint

  # Only genuinely completed scientific outcomes are terminal. Policy-stopped
  # and interrupted rows are reset to pending so a later resume can continue
  # the unfinished study.
  terminal_mask <- is_terminal_task_status(checkpoint_grid$status)

  if (!any(terminal_mask)) {
    # No terminal tasks to merge
    return(fresh_grid)
  }

  terminal_tasks <- checkpoint_grid[terminal_mask, ]

  # Build lookup for fast status retrieval
  status_lookup <- stats::setNames(
    terminal_tasks$status,
    terminal_tasks$task_id
  )

  # Update status for matching task IDs
  # Use a single match() pass over the lookup (O(n + m)) instead of an
  # %in% scan per row (O(n * m)) — large task grids make the difference
  # quadratic.
  task_ids <- fresh_grid$task_id
  hit <- match(task_ids, names(status_lookup))
  hit_ids <- which(!is.na(hit))
  fresh_grid$status[hit_ids] <- unname(status_lookup)[hit[hit_ids]]

  fresh_grid
}

#' Merge Prior and New Results
#'
#' Combines results from a resumed run with new results from continued
#' execution. Duplicate task rows must be identical or an integrity error
#' is raised.
#'
#' @param prior_results Data frame of results from checkpoint.
#' @param new_results Data frame of results from new execution.
#'
#' @return Combined data frame with unique task_id values.
#'
#' @details
#' If both inputs are NULL or empty, returns NULL. If only one input
#' has data, returns that input. Otherwise, removes any task_ids from
#' prior_results that appear in new_results, then combines.
#'
#' @keywords internal
merge_results <- function(prior_results, new_results) {
  # Handle empty/NULL cases
  if (is.null(prior_results) || nrow(prior_results) == 0) {
    return(new_results)
  }

  if (is.null(new_results) || nrow(new_results) == 0) {
    return(prior_results)
  }

  duplicate_ids <- intersect(prior_results$task_id, new_results$task_id)

  if (length(duplicate_ids) > 0) {
    for (task_id in duplicate_ids) {
      prior_row <- prior_results[
        prior_results$task_id == task_id,
        ,
        drop = FALSE
      ]
      new_row <- new_results[new_results$task_id == task_id, , drop = FALSE]

      if (nrow(prior_row) != 1 || nrow(new_row) != 1) {
        cli::cli_abort(
          "Duplicate terminal rows detected for task_id '{task_id}'",
          class = "bayesim_checkpoint_error"
        )
      }

      # Policy-stopped/interrupted rows are placeholders, not terminal
      # outcomes. A resumed execution is expected to replace them with the
      # newly computed result.
      prior_status <- as.character(prior_row$status[[1]])
      if (is_resumable_task_status(prior_status)) {
        next
      }

      if (!rows_equivalent(prior_row, new_row)) {
        cli::cli_abort(
          paste(
            "Conflicting duplicate terminal rows detected for task_id",
            shQuote(task_id)
          ),
          class = "bayesim_checkpoint_error"
        )
      }
    }
  }

  prior_only <- prior_results[
    !prior_results$task_id %in% duplicate_ids,
    ,
    drop = FALSE
  ]

  # Checkpoints from adjacent schema versions may differ only by an all-NA
  # optional column. Align the union before binding so an equivalent duplicate
  # can be replaced by the new row without making unrelated prior rows
  # impossible to retain. Assign row-count-length NA vectors: a bare NA errors
  # when the receiving side is a 0-row frame, which happens when the resumed
  # execution re-covers every prior task (resume-to-completion; #63).
  all_cols <- union(names(prior_only), names(new_results))
  for (col in setdiff(all_cols, names(prior_only))) {
    prior_only[[col]] <- rep(NA, nrow(prior_only))
  }
  for (col in setdiff(all_cols, names(new_results))) {
    new_results[[col]] <- rep(NA, nrow(new_results))
  }
  combined <- rbind(
    prior_only[, all_cols, drop = FALSE],
    new_results[, all_cols, drop = FALSE]
  )

  # Reset row names for clean output
  rownames(combined) <- NULL

  combined
}

normalize_result_row <- function(x) {
  # S1: timing is never part of the task's identity — a re-run of an identical
  # task legitimately produces a different wall-clock time. Excluding it keeps
  # the resume duplicate-check from aborting on spurious timing differences.
  x <- x[, setdiff(names(x), "timing_total"), drop = FALSE]
  cols <- sort(names(x))
  x <- x[, cols, drop = FALSE]

  lapply(x, function(col) {
    if (is.factor(col)) {
      as.character(col)
    } else {
      col
    }
  })
}

# Compare two result rows for equivalence, tolerating benign representation
# differences that `identical()` would flag spuriously:
# - NA type mismatches (e.g. logical NA vs NA_real_): results_to_dataframe()
#   fills missing columns with logical NA, so the same all-missing value can
#   legitimately be represented with different NA types across checkpoints.
# - Column sets that differ only by NA-filled columns: a column absent from
#   one row is treated as all-NA.
# Rows agree when every column over the union of names has an identical NA
# pattern and identical non-NA values. All-NA values of any type are equal.
rows_equivalent <- function(x, y) {
  x <- normalize_result_row(x)
  y <- normalize_result_row(y)

  cols <- union(names(x), names(y))
  for (col in cols) {
    vx <- x[[col]]
    vy <- y[[col]]
    if (is.null(vx)) {
      vx <- NA
    }
    if (is.null(vy)) {
      vy <- NA
    }

    vx_na <- is.na(vx)
    vy_na <- is.na(vy)
    if (!identical(vx_na, vy_na)) {
      return(FALSE)
    }
    if (!all(vx_na) && !identical(vx[!vx_na], vy[!vy_na])) {
      return(FALSE)
    }
  }
  TRUE
}

rehydrate_config_from_manifest <- function(result_path) {
  manifest <- read_run_manifest(result_path)

  if (is.null(manifest) || is.null(manifest$config_spec)) {
    cli::cli_abort(
      "Run manifest does not contain a rehydratable configuration; please supply config explicitly"
    )
  }

  spec <- manifest$config_spec
  # Runtime policy is mutable across compatible resumes. The latest valid
  # checkpoint is its atomic source of truth; the manifest remains a legacy
  # fallback and the immutable home of scientific StudySpec inputs.
  checkpoint <- get_latest_valid_checkpoint(
    result_path,
    config_fingerprint = manifest$config_fingerprint,
    load_outcomes = FALSE
  )
  policy <- checkpoint$meta$run_policy %||% manifest$run_policy %||% list()
  data_generator <- rehydrate_function_spec(spec$data_generator_spec)
  fitter <- rehydrate_s7_spec(spec$fitter_spec)
  metrics <- lapply(spec$metrics_spec %||% list(), rehydrate_s7_spec)

  data_grid <- normalize_manifest_df(spec$data_grid)
  fit_grid <- normalize_manifest_df(spec$fit_grid)
  task_grid <- normalize_manifest_df(spec$task_grid)
  retain <- normalize_manifest_retain(
    policy$retain %||% manifest$retention_spec %||% spec$retain
  )
  max_errors <- normalize_manifest_numeric(
    policy$max_errors %||% spec$max_errors,
    default = Inf
  )
  stop_on <- normalize_manifest_stop_on(policy$stop_on %||% spec$stop_on)
  checkpoint_every <- as.integer(normalize_manifest_numeric(
    policy$checkpoint_every %||% spec$checkpoint_every,
    default = 50L
  ))
  keep_checkpoints <- as.integer(normalize_manifest_numeric(
    policy$keep_checkpoints %||% spec$keep_checkpoints,
    default = 2L
  ))
  seed <- as.integer(normalize_manifest_numeric(spec$seed, default = 1L))
  n_replicates <- as.integer(normalize_manifest_numeric(
    spec$n_replicates,
    default = 1L
  ))

  simulation_config(
    data_grid = data_grid,
    fit_grid = fit_grid,
    task_grid = task_grid,
    data_generator = data_generator,
    fitter = fitter,
    metrics = metrics,
    n_replicates = n_replicates,
    seed = seed,
    result_path = result_path,
    checkpoint_format = policy$checkpoint_format %||%
      manifest$checkpoint_format %||%
      spec$checkpoint_format %||%
      "rds",
    checkpoint_every = checkpoint_every,
    keep_checkpoints = keep_checkpoints,
    retain = retain,
    max_errors = max_errors,
    stop_on = stop_on
  )
}

# F6: Truthful resume guidance --------------------------------------------
#
# print.bayesim_simulation_result() used to advertise the configless
# `resume_simulation("<path>")` command unconditionally, but that command only
# works when every executable component recorded in the run manifest (data
# generator, fitter, metrics) is namespaced and therefore rehydratable. A
# normal exported fixed_truth_generator() returns a closure whose manifest
# executable is non-rehydratable, so the advertised command would be rejected
# by rehydrate_function_spec()/rehydrate_s7_spec(). Inspect the actual
# manifest before recommending anything.

# Assess whether a run's manifest supports configless resume.
#
# Returns a list with `rehydratable` (logical) and `components` (character
# vector naming the executable components the manifest cannot restore).
# An unreadable manifest, or one without a config_spec, is reported as
# non-rehydratable: configless resume would fail for it too, so truthful
# guidance must ask for the original config.
run_manifest_rehydration_status <- function(result_path) {
  status <- list(rehydratable = FALSE, components = character())

  manifest <- tryCatch(
    read_run_manifest(result_path),
    error = function(e) NULL
  )
  if (is.null(manifest) || is.null(manifest$config_spec)) {
    return(status)
  }

  spec <- manifest$config_spec
  component_rehydratable <- function(entry) {
    is.list(entry) && length(entry) > 0L && isTRUE(entry$rehydratable)
  }

  components <- character()
  if (!component_rehydratable(spec$data_generator_spec)) {
    components <- c(components, "data_generator")
  }
  if (!component_rehydratable(spec$fitter_spec)) {
    components <- c(components, "fitter")
  }
  metric_specs <- if (is.list(spec$metrics_spec)) spec$metrics_spec else list()
  metric_ok <- vapply(metric_specs, component_rehydratable, logical(1))
  if (!all(metric_ok)) {
    # Name the offending classes when available so the user knows which
    # metric objects must be rebuilt.
    bad <- metric_specs[!metric_ok]
    classes <- vapply(
      bad,
      function(entry) {
        cl <- entry$class %||% NA_character_
        if (isTRUE(is.na(cl))) "unknown class" else cl
      },
      character(1)
    )
    components <- c(
      components,
      paste0("metrics (", paste(unique(classes), collapse = ", "), ")")
    )
  }

  status$components <- components
  status$rehydratable <- length(components) == 0L
  status
}

# Build the resume guidance lines printed for a completed run.
#
# Only recommend the configless `resume_simulation(path)` command when the
# manifest inspection says it will work. Otherwise print an equally exact
# command that requires the original config, plus the reason (the manifest
# cannot rehydrate script-defined closures or non-namespaced objects, and
# resume does not pretend otherwise).
resume_guidance_lines <- function(result_path) {
  status <- run_manifest_rehydration_status(result_path)

  if (isTRUE(status$rehydratable)) {
    return(paste0(
      "  Resume with: resume_simulation(",
      encodeString(result_path, quote = "\""),
      ")"
    ))
  }

  lines <- paste0(
    "  Resume with: resume_simulation(",
    encodeString(result_path, quote = "\""),
    ", config = config)"
  )
  hint <- paste0(
    "  (`config` is the original SimulationConfig object ",
    "used to create this run)"
  )
  reason <- if (length(status$components) > 0L) {
    paste0(
      "  Manifest cannot restore: ",
      paste(status$components, collapse = ", "),
      "; supply the config used to create this run."
    )
  } else {
    paste0(
      "  No rehydratable configuration found in the manifest; ",
      "supply the config used to create this run."
    )
  }
  c(lines, hint, reason)
}

rehydrate_function_spec <- function(spec) {
  if (is.null(spec) || !isTRUE(spec$rehydratable) || is.null(spec$reference)) {
    cli::cli_abort(
      "Run manifest references non-rehydratable executable components; please supply config explicitly"
    )
  }

  validate_namespace_version(spec$reference$package, spec$reference$version)

  get(spec$reference$name, envir = asNamespace(spec$reference$package))
}

rehydrate_s7_spec <- function(spec) {
  if (is.null(spec) || (length(spec) == 1 && is.na(spec))) {
    return(NULL)
  }

  if (!isTRUE(spec$rehydratable) || is.null(spec$class)) {
    cli::cli_abort(
      "Run manifest references non-rehydratable executable components; please supply config explicitly"
    )
  }

  class_parts <- strsplit(spec$class, "::", fixed = TRUE)[[1]]
  if (length(class_parts) != 2) {
    cli::cli_abort(
      "Cannot rehydrate S7 object without namespaced class reference"
    )
  }

  validate_namespace_version(class_parts[[1]], spec$package_version)

  constructor <- get(class_parts[[2]], envir = asNamespace(class_parts[[1]]))
  props <- spec$properties %||% list()
  props <- Map(normalize_manifest_property, props, names(props))
  do.call(constructor, props)
}

normalize_manifest_property <- function(value, name = NULL) {
  if (identical(name, "needs") && is.list(value) && length(value) == 0L) {
    return(character())
  }
  if (
    is.list(value) &&
      is.null(names(value)) &&
      all(vapply(value, function(x) !is.list(x) && length(x) == 1L, logical(1)))
  ) {
    return(unlist(value, use.names = FALSE))
  }
  value
}

validate_namespace_version <- function(package, expected_version = NULL) {
  if (is.null(expected_version)) {
    return(invisible(TRUE))
  }

  current_version <- as.character(getNamespaceVersion(package))
  if (!identical(current_version, expected_version)) {
    cli::cli_abort(
      c(
        "Component version mismatch during resume rehydration",
        "x" = "Package '{package}' version is '{current_version}', expected '{expected_version}'"
      )
    )
  }

  invisible(TRUE)
}

normalize_manifest_df <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.list(x) && length(x) == 0) {
    return(NULL)
  }

  if (is.data.frame(x)) {
    return(x)
  }

  if (is.list(x) && all(vapply(x, is.list, logical(1)))) {
    # Combine list-of-lists into a data.frame, filling missing columns with NA
    # (replaces dplyr::bind_rows to drop the dplyr dependency).
    all_cols <- unique(unlist(lapply(x, names)))
    rows <- lapply(x, function(r) {
      missing <- setdiff(all_cols, names(r))
      if (length(missing) > 0L) {
        r[missing] <- rep(NA, length(missing))
      }
      as.data.frame(r[all_cols], stringsAsFactors = FALSE)
    })
    return(do.call(rbind, rows))
  }

  tryCatch(as.data.frame(x), error = function(e) x)
}

normalize_manifest_retain <- function(x) {
  if (is.null(x) || (is.list(x) && length(x) == 0)) {
    return(c("metrics", "diagnostics"))
  }

  if (is.character(x)) {
    return(x)
  }

  if (!is.list(x)) {
    return(c("metrics", "diagnostics"))
  }

  lapply(x, function(value) {
    if (is.character(value)) {
      value
    } else if (is.list(value)) {
      unlist(value, use.names = FALSE)
    } else {
      as.character(value)
    }
  })
}

normalize_manifest_numeric <- function(x, default = NULL) {
  if (is.null(x) || (is.list(x) && length(x) == 0)) {
    return(default)
  }

  if (is.character(x) && identical(x, "Inf")) {
    return(Inf)
  }

  as.numeric(x)
}

normalize_manifest_stop_on <- function(x) {
  if (is.null(x) || (is.list(x) && length(x) == 0L)) {
    return(NULL)
  }
  if (!is.list(x)) {
    return(x)
  }
  Map(normalize_manifest_property, x, names(x))
}
