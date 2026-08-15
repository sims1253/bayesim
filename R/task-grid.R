#' @keywords internal
#' @importFrom parallel nextRNGStream
NULL

#' Create Task RNG Streams
#'
#' Precomputes deterministic RNG streams for each task using L'Ecuyer-CMRG.
#' Each task receives its own independent RNG stream, ensuring reproducibility
#' regardless of execution order (sequential, parallel, or resumed).
#'
#' @param global_seed Integer. The base seed for the simulation.
#' @param n_tasks Integer. Number of tasks to generate streams for.
#'
#' @return A list of length `n_tasks`, where each element is an integer vector
#'   representing the `.Random.seed` state for that task.
#'
#' @keywords internal
#'
#' @details The caller's RNG kind and seed are restored on exit, including when
#'   `.Random.seed` did not exist before the call.
create_task_rng_streams <- function(global_seed, n_tasks) {
  old_kind <- base::RNGkind()
  seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (seed_existed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit(
    {
      do.call(base::RNGkind, as.list(old_kind))
      if (seed_existed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (
        exists(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE
        )
      ) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )

  base::RNGkind("L'Ecuyer-CMRG")
  set.seed(global_seed)

  streams <- vector("list", n_tasks)
  seed <- .Random.seed
  for (i in seq_len(n_tasks)) {
    streams[[i]] <- seed
    seed <- parallel::nextRNGStream(seed)
  }
  streams
}

task_id_widths <- function(data_idx, fit_idx, rep_idx) {
  list(
    data = max(3L, nchar(as.character(max(data_idx, na.rm = TRUE)))),
    fit = max(3L, nchar(as.character(max(fit_idx, na.rm = TRUE)))),
    rep = max(5L, nchar(as.character(max(rep_idx, na.rm = TRUE))))
  )
}

#' Create task ID from indices
#'
#' @param data_idx Integer. Data index.
#' @param fit_idx Integer. Fit index.
#' @param rep_idx Integer. Replication index.
#' @param widths Optional list with data, fit, rep widths. Auto-computed if NULL.
#'
#' @return Character string task ID.
#' @keywords internal
make_task_id <- function(data_idx, fit_idx, rep_idx, widths = NULL) {
  if (is.null(widths)) {
    widths <- task_id_widths(data_idx, fit_idx, rep_idx)
  }

  sprintf(
    paste0(
      "d%0",
      widths$data,
      "d_f%0",
      widths$fit,
      "d_r%0",
      widths$rep,
      "d"
    ),
    data_idx,
    fit_idx,
    rep_idx
  )
}

canonicalize_task_grid <- function(task_grid, config) {
  grid <- tibble::as_tibble(task_grid)

  if ("task_id" %in% names(grid) && anyDuplicated(grid$task_id)) {
    cli::cli_abort("task_grid contains duplicate task_id values")
  }

  if (!"rep_idx" %in% names(grid)) {
    grid$rep_idx <- rep.int(1L, nrow(grid))
  }
  grid$rep_idx <- as.integer(grid$rep_idx)

  has_explicit_specs <- all(c("data_spec", "fit_spec") %in% names(grid))

  if (has_explicit_specs) {
    if (!is.list(grid$data_spec) || !is.list(grid$fit_spec)) {
      cli::cli_abort(
        "task_grid$data_spec and task_grid$fit_spec must be list-columns"
      )
    }

    if (!"data_idx" %in% names(grid)) {
      grid$data_idx <- seq_len(nrow(grid))
    }
    if (!"fit_idx" %in% names(grid)) {
      grid$fit_idx <- rep.int(1L, nrow(grid))
    }
  } else {
    if (!all(c("data_idx", "fit_idx") %in% names(grid))) {
      cli::cli_abort("task_grid must contain data_idx and fit_idx columns")
    }

    grid$data_idx <- as.integer(grid$data_idx)
    grid$fit_idx <- as.integer(grid$fit_idx)

    if (anyNA(grid$data_idx) || any(grid$data_idx < 1L)) {
      cli::cli_abort("task_grid$data_idx must contain positive integers")
    }
    if (anyNA(grid$fit_idx) || any(grid$fit_idx < 1L)) {
      cli::cli_abort("task_grid$fit_idx must contain positive integers")
    }

    if (is.null(config@data_grid) || is.null(config@fit_grid)) {
      cli::cli_abort(
        paste(
          "task_grid with data_idx/fit_idx requires data_grid and fit_grid in the configuration"
        )
      )
    }
  }

  grid <- grid[order(grid$data_idx, grid$fit_idx, grid$rep_idx), , drop = FALSE]
  rownames(grid) <- NULL

  widths <- task_id_widths(grid$data_idx, grid$fit_idx, grid$rep_idx)
  grid$task_id <- make_task_id(
    grid$data_idx,
    grid$fit_idx,
    grid$rep_idx,
    widths
  )
  if (anyDuplicated(grid$task_id)) {
    cli::cli_abort(
      "task_grid contains duplicate task identities; data_idx, fit_idx, and rep_idx must be unique"
    )
  }
  grid$rng_seed <- I(create_task_rng_streams(config@seed, nrow(grid)))
  grid$status <- "pending"
  grid$stop_reason <- NA_character_

  tibble::as_tibble(grid)
}

#' Create Task Grid from Configuration
#'
#' Generates the complete, deterministic task table for a simulation.
#' Tasks are ordered lexicographically by task_id for reproducibility.
#'
#' @param config A SimulationConfig S7 object.
#'
#' @return A tibble with one row per task, containing columns:
#'   - `task_id`: Character identifier in format "dXXX_fXXX_rXXXXX"
#'   - `data_idx`: Integer index into the data_grid (1-based)
#'   - `fit_idx`: Integer index into the fit_grid (1-based)
#'   - `rep_idx`: Integer replicate number (1 to n_replicates)
#'   - `rng_seed`: List column containing precomputed RNG stream for each task
#'   - `status`: Character status, initialized to "pending"
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' config <- simulation_config(
#'   data_grid = data.frame(n = c(100, 500)),
#'   fit_grid = data.frame(model = c("baseline", "full")),
#'   data_generator = my_data_gen,
#'   n_replicates = 10L,
#'   seed = 42L
#' )
#'
#' task_grid <- create_task_grid(config)
#' # task_grid has 2 * 2 * 10 = 40 rows
#' }
create_task_grid <- function(config) {
  # Validate config
  if (!is_simulation_config(config)) {
    cli::cli_abort("config must be a SimulationConfig object")
  }

  if (!is.null(config@task_grid)) {
    return(canonicalize_task_grid(config@task_grid, config))
  }

  n_data <- nrow(config@data_grid)
  n_fit <- nrow(config@fit_grid)
  n_rep <- config@n_replicates

  total_tasks <- n_data * n_fit * n_rep

  # Generate all combinations
  grid <- expand.grid(
    data_idx = seq_len(n_data),
    fit_idx = seq_len(n_fit),
    rep_idx = seq_len(n_rep),
    KEEP.OUT.ATTRS = FALSE
  )

  # Sort lexicographically for deterministic ordering
  grid <- grid[order(grid$data_idx, grid$fit_idx, grid$rep_idx), ]

  # Reset row names after sorting
  rownames(grid) <- NULL

  widths <- task_id_widths(grid$data_idx, grid$fit_idx, grid$rep_idx)
  grid$task_id <- make_task_id(
    grid$data_idx,
    grid$fit_idx,
    grid$rep_idx,
    widths
  )

  # Precompute RNG streams
  rng_streams <- create_task_rng_streams(config@seed, total_tasks)
  grid$rng_seed <- I(rng_streams)

  # Initialize status
  grid$status <- "pending"
  grid$stop_reason <- NA_character_

  tibble::as_tibble(grid)
}

#' Get a task specification by precomputed row index
#'
#' The execution loop resolves a whole batch of task IDs once and calls this
#' helper by index, avoiding repeated full-grid scans for large studies.
#'
#' @param task_grid A task grid tibble.
#' @param row_idx Scalar row index.
#' @param config The SimulationConfig used to create the grid.
#' @return A named task specification list containing the task and grid
#'   indices, data and fit specifications, task context, and RNG seed.
#' @keywords internal
get_task_spec_at <- function(task_grid, row_idx, config) {
  row <- task_grid[row_idx, , drop = FALSE]
  task_id <- row$task_id[[1]]

  list(
    task_id = task_id,
    data_idx = row$data_idx[[1]],
    fit_idx = row$fit_idx[[1]],
    rep_idx = row$rep_idx[[1]],
    data_spec = if ("data_spec" %in% names(row)) {
      row$data_spec[[1]]
    } else {
      as.list(config@data_grid[row$data_idx[[1]], , drop = FALSE])
    },
    fit_spec = if ("fit_spec" %in% names(row)) {
      row$fit_spec[[1]]
    } else {
      as.list(config@fit_grid[row$fit_idx[[1]], , drop = FALSE])
    },
    task_ctx = list(
      task_id = task_id,
      data_idx = row$data_idx[[1]],
      fit_idx = row$fit_idx[[1]],
      rep_idx = row$rep_idx[[1]]
    ),
    rng_seed = row$rng_seed[[1]]
  )
}

#' Filter Tasks by Status
#'
#' Filters a task grid to include only tasks with specified statuses.
#'
#' @param task_grid A task grid tibble.
#' @param status Character vector of statuses to include.
#'   Valid statuses: "pending", "success", "failed", "skipped", or
#'   "cancelled".
#'
#' @return A filtered task grid tibble.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Get all completed or skipped tasks
#' done <- filter_tasks_by_status(task_grid, c("success", "skipped"))
#' }
filter_tasks_by_status <- function(task_grid, status) {
  task_grid[task_grid$status %in% status, ]
}

#' Get Pending Tasks
#'
#' Convenience function to retrieve all tasks that have not yet been processed.
#'
#' @param task_grid A task grid tibble.
#'
#' @return A task grid tibble containing only tasks with status "pending".
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' pending <- get_pending_tasks(task_grid)
#' nrow(pending)  # Number of remaining tasks
#' }
get_pending_tasks <- function(task_grid) {
  filter_tasks_by_status(task_grid, "pending")
}

#' Validate Task ID Format
#'
#' Checks if a string is a valid task ID in the format "dXXX_fXXX_rXXXXX".
#'
#' @param task_id Character string to validate.
#'
#' @return TRUE if valid, FALSE otherwise.
#'
#' @keywords internal
validate_task_id <- function(task_id) {
  grepl("^d[0-9]+_f[0-9]+_r[0-9]+$", task_id)
}

#' Parse Task ID Components
#'
#' Extracts the data_idx, fit_idx, and rep_idx from a task ID string.
#'
#' @param task_id Character task ID in format "dXXX_fXXX_rXXXXX".
#'
#' @return Named integer vector with elements data_idx, fit_idx, rep_idx,
#'   or NULL if the task_id format is invalid.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' components <- parse_task_id("d002_f003_r00042")
#' # components$data_idx == 2
#' # components$fit_idx == 3
#' # components$rep_idx == 42
#' }
parse_task_id <- function(task_id) {
  if (!validate_task_id(task_id)) {
    return(NULL)
  }

  parts <- strsplit(task_id, "_")[[1]]

  list(
    data_idx = as.integer(sub("^d", "", parts[1])),
    fit_idx = as.integer(sub("^f", "", parts[2])),
    rep_idx = as.integer(sub("^r", "", parts[3]))
  )
}

#' Get Task Count Summary
#'
#' Returns a summary of task counts by status.
#'
#' @param task_grid A task grid tibble.
#'
#' @return Named integer vector with counts for each status.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' summary <- get_task_count_summary(task_grid)
#' # summary["pending"]  # Number of pending tasks
#' # summary["success"]  # Number of successful tasks
#' }
get_task_count_summary <- function(task_grid) {
  # Preserve the compact historical shape for ordinary ledgers while exposing
  # explicit lifecycle states whenever they actually occur.
  statuses <- c("pending", TASK_TERMINAL_STATUSES, "skipped")
  statuses <- c(
    statuses,
    intersect(
      setdiff(TASK_RESUMABLE_STATUSES, statuses),
      unique(task_grid$status)
    )
  )
  counts <- table(factor(task_grid$status, levels = statuses))
  stats::setNames(as.integer(counts), statuses)
}
