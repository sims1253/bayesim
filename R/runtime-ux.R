# Runtime UX helpers (Workstream F): preflight, failure surfacing,
# print/as_tibble polish. These are additive, exported utilities.

# F1: preflight -----------------------------------------------------------

#' Preflight check for a simulation configuration
#'
#' @description
#' Inspects a configuration *before* a long run and reports: total task count,
#' grid dimensions, metrics and their `needs` vs the fitter's capabilities
#' (surfacing mismatch warnings before a 10-hour run), whether mirai daemons
#' are set, and (for BrmsFitter) the number of distinct models to compile.
#'
#' When run automatically at the top of [run_simulation()], it prints a
#' condensed one-line summary.
#'
#' @param config A `SimulationConfig`.
#' @param pilot Logical; if TRUE, run a single pilot task and extrapolate its
#'   wall-clock time to an estimate of the total study time (default FALSE). The
#'   estimate is printed and returned in the `pilot_seconds` and
#'   `estimated_total_seconds` elements of the returned list. A pilot that
#'   cannot execute (e.g. a Stan fitter with no CmdStan installed) yields NA
#'   estimates rather than failing the preflight. Note that for Stan fitters a
#'   pilot compiles the model, so it is not instant.
#' @param condensed Logical; if TRUE, print a single one-line summary instead of
#'   the full report (used by run_simulation).
#' @return Invisibly, a named list of preflight information.
#' @export
#' @seealso [run_simulation()], [simulation_config()]
#' @examples
#' \dontrun{
#' preflight(config)
#' }
preflight <- function(config, pilot = FALSE, condensed = FALSE) {
  if (!is_simulation_config(config)) {
    stop(bayesim_config_error("preflight requires a SimulationConfig"))
  }
  n_tasks <- get_total_tasks(config)
  n_data <- if (!is.null(config@data_grid)) nrow(config@data_grid) else NA
  n_fit <- if (!is.null(config@fit_grid)) nrow(config@fit_grid) else NA
  n_rep <- config@n_replicates
  metrics <- config@metrics %||% list()
  metric_names <- vapply(metrics, function(m) m@name, character(1))
  metric_needs <- unique(unlist(lapply(metrics, function(m) m@needs)))

  fitter <- config@fitter
  caps <- if (!is.null(fitter) && S7::S7_inherits(fitter)) {
    c(
      if (isTRUE(fitter@supports_predictions)) "predictions",
      if (isTRUE(fitter@supports_log_lik)) "log_lik",
      if (isTRUE(fitter@supports_loo)) "loo",
      if (isTRUE(fitter@supports_epred)) "epred"
    )
  } else {
    character()
  }
  unmet <- setdiff(metric_needs, caps)

  daemons_set <- isTRUE(mirai::daemons_set())
  n_compile <- if (
    S7::S7_inherits(fitter, BrmsFitter) && !is.null(config@fit_grid)
  ) {
    hashes <- vapply(
      seq_len(nrow(config@fit_grid)),
      function(i) {
        spec <- model_spec_from_grid_row(config@fit_grid, i)
        model_spec_hash(
          spec$formula,
          spec$family,
          spec$prior,
          spec$stanvars,
          fitter@backend
        )
      },
      character(1)
    )
    length(unique(hashes))
  } else if (S7::S7_inherits(fitter, CmdStanFitter_class)) {
    1L
  } else {
    0L
  }

  info <- list(
    n_tasks = n_tasks,
    data_grid = n_data,
    fit_grid = n_fit,
    n_replicates = n_rep,
    metrics = metric_names,
    metric_needs = metric_needs,
    fitter_capabilities = caps,
    unmet_needs = unmet,
    daemons_set = daemons_set,
    n_compile = n_compile,
    pilot_seconds = NA_real_,
    estimated_total_seconds = NA_real_
  )

  # R4: opt-in pilot timing. Run the first grid cell once and extrapolate.
  if (pilot) {
    pilot_seconds <- run_pilot_task(config)
    if (!is.null(pilot_seconds)) {
      info$pilot_seconds <- pilot_seconds
      info$estimated_total_seconds <- pilot_seconds * n_tasks
    }
  }

  if (condensed) {
    # mirai::daemons() (no args) is not an accessor in mirai >= 2.x — it errors
    # with "argument n is missing". Use status()$connections for the count.
    workers_str <- if (daemons_set) {
      n_daemons <- tryCatch(mirai::status()$connections, error = function(e) {
        NULL
      })
      if (length(n_daemons)) paste0("; ", n_daemons, " daemons") else ""
    } else {
      ""
    }
    compile_str <- if (n_compile > 0L) {
      paste0("; ", n_compile, " models to compile")
    } else {
      ""
    }
    cli::cli_inform(
      "{n_tasks} tasks = {n_data} data x {n_fit} fit x {n_rep} reps{compile_str}{workers_str}"
    )
    # R1a: surface unmet metric needs in condensed mode too — this is the path
    # run_simulation() uses by default, so an all-NA study is never a silent
    # surprise.
    if (length(unmet)) {
      cli::cli_warn(c(
        "Metrics need capabilities the fitter does not support.",
        i = paste("unmet needs:", paste(unmet, collapse = ", "))
      ))
    }
    if (pilot) print_pilot_estimate(info)
  } else {
    cli::cli_h1("bayesim preflight")
    cli::cli_text(
      "{n_tasks} tasks = {n_data} data x {n_fit} fit x {n_rep} replicates"
    )
    if (n_compile > 0L) {
      cli::cli_text("Models to compile: {n_compile}")
    }
    cli::cli_text("Metrics: {.val {metric_names}}")
    if (length(unmet)) {
      cli::cli_warn(c(
        "Metrics need capabilities the fitter does not support.",
        i = paste("unmet needs:", paste(unmet, collapse = ", "))
      ))
    } else if (length(metric_needs)) {
      cli::cli_text(
        "Metric needs ({.val {metric_needs}}) all satisfied by the fitter."
      )
    }
    cli::cli_text("Daemons set: {daemons_set}")
    if (pilot) print_pilot_estimate(info)
  }

  invisible(info)
}

# Print the pilot timing estimate (R4). Shared by the condensed and full
# preflight reports.
print_pilot_estimate <- function(info) {
  if (is.finite(info$estimated_total_seconds)) {
    cli::cli_text(
      "Pilot task: {round(info$pilot_seconds, 3)}s; estimated total ~{round(info$estimated_total_seconds, 1)}s"
    )
  } else {
    cli::cli_text("(pilot timing unavailable for this configuration)")
  }
  invisible(NULL)
}

# R4: run a single representative task (the first grid cell) to measure
# per-task wall-clock time. Returns seconds (numeric scalar) or NULL when the
# pilot cannot execute for this configuration (e.g. a Stan fitter with no
# CmdStan installed). Used by preflight(pilot = TRUE). Runs sequentially on the
# controller with no result_path, so no metric artifacts leak into a real
# results directory, and restores the caller's RNG state on exit.
run_pilot_task <- function(config) {
  if (!is.null(config@task_grid)) {
    # The custom task_grid may carry data_spec/fit_spec list-columns, plain
    # data_idx/fit_idx indices, or (after canonicalization) both. Handle each.
    if ("data_spec" %in% names(config@task_grid)) {
      data_spec <- config@task_grid$data_spec[[1L]]
      data_idx <- config@task_grid$data_idx[[1L]] %||% 1L
    } else if (!is.null(config@data_grid)) {
      data_idx <- config@task_grid$data_idx[[1L]] %||% 1L
      data_spec <- as.list(config@data_grid[data_idx, , drop = FALSE])
    } else {
      return(NULL)
    }
    if ("fit_spec" %in% names(config@task_grid)) {
      fit_spec <- config@task_grid$fit_spec[[1L]]
      fit_idx <- config@task_grid$fit_idx[[1L]] %||% 1L
    } else if (!is.null(config@fit_grid)) {
      fit_idx <- config@task_grid$fit_idx[[1L]] %||% 1L
      fit_spec <- as.list(config@fit_grid[fit_idx, , drop = FALSE])
    } else {
      return(NULL)
    }
  } else {
    data_idx <- 1L
    fit_idx <- 1L
    data_spec <- as.list(config@data_grid[1L, , drop = FALSE])
    fit_spec <- as.list(config@fit_grid[1L, , drop = FALSE])
  }

  # Restore the caller's RNG after the pilot (run_task_safe advances the global
  # .Random.seed via the task stream).
  old_seed <- if (
    exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  ) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit(
    {
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    },
    add = TRUE
  )

  task_ctx <- list(
    task_id = "d001_f001_r00001",
    data_idx = data_idx,
    fit_idx = fit_idx,
    rep_idx = 1L
  )
  streams <- create_task_rng_streams(config@seed, 1L)

  task <- list(
    task_id = task_ctx$task_id,
    data_spec = data_spec,
    fit_spec = fit_spec,
    task_ctx = task_ctx,
    rng_seed = streams[[1L]]
  )

  config_spec <- as_config_spec(config)
  config_spec$data_generator <- config@data_generator
  config_spec$result_path <- NULL
  config_spec$package_name <- utils::packageName()

  result <- tryCatch(
    run_task_safe(
      task,
      config_spec,
      config@fitter,
      config@metrics,
      retain = config@retain
    ),
    error = function(e) NULL
  )
  if (is.null(result) || !identical(result$status, "success")) {
    return(NULL)
  }
  result$timing$total
}

# F2: failure surfacing ---------------------------------------------------

#' Extract failed tasks from a simulation result
#'
#' @description Returns the failed-task errors joined with grid columns, as a
#' tibble — convenient for diagnosing what went wrong after a run with
#' failures. Also used by run_simulation to print a compact failure summary.
#' @param result A `bayesim_simulation_result`.
#' @return A tibble with grid columns plus `error_class`, `error_message`.
#' @export
#' @seealso [run_simulation()]
#' @examples
#' \dontrun{
#' failed_tasks(result)
#' }
failed_tasks <- function(result) {
  df <- if (is.data.frame(result)) result else result$summary
  if (is.null(df) || !is.data.frame(df)) {
    stop(bayesim_config_error(
      "failed_tasks requires a simulation result or summary"
    ))
  }
  failed <- df[!is.na(df$status) & df$status == "failed", , drop = FALSE]
  if (nrow(failed) == 0L) {
    return(tibble::tibble())
  }
  cols <- intersect(
    c(
      "task_id",
      "error_class",
      "error_message",
      grep("^(data_|fit_|rep_idx)", names(failed), value = TRUE)
    ),
    names(failed)
  )
  tibble::as_tibble(failed[, cols, drop = FALSE])
}

# Compact failure summary printer (used internally by run_simulation).
print_failure_summary <- function(result) {
  ft <- failed_tasks(result)
  if (nrow(ft) == 0L) {
    return(invisible(NULL))
  }
  cli::cli_alert_warning("{nrow(ft)} task(s) failed:")
  by_class <- split(ft, ft$error_class)
  for (cl in names(by_class)) {
    sub <- by_class[[cl]]
    first_msg <- sub$error_message[1]
    extra <- if (nrow(sub) > 1L) {
      paste0(" (and ", nrow(sub) - 1L, " more)")
    } else {
      ""
    }
    cli::cli_text("  {.val {cl}}: {first_msg}{extra}")
  }
  invisible(NULL)
}

# F7: end-of-run summary ----------------------------------------------------

# Human phrasing for the stop_reason values execute_tasks() records on tasks
# that were never executed (policy-stopped work is resumable by design).
format_stop_reasons <- function(reasons) {
  labels <- c(
    max_errors = "error budget exhausted",
    adaptive_stop = "adaptive stopping target reached"
  )
  known <- labels[reasons]
  unknown <- is.na(known)
  if (any(unknown)) {
    known[unknown] <- reasons[unknown]
  }
  paste(unique(known), collapse = "; ")
}

# Print the single end-of-run block: completion status, task counts, stop
# reason, results path, and — when unexecuted work remains — the literal
# resume command (truthful about configless resume, see resume_guidance_lines),
# followed by the per-class failure detail. Called by run_simulation() when
# verbose = TRUE; independent of the task progress bar.
print_run_summary <- function(result) {
  grid <- result$task_grid
  statuses <- if (is.data.frame(grid) && nrow(grid) > 0L) {
    grid$status
  } else {
    vapply(
      result$task_results,
      function(tr) if (is.null(tr)) "pending" else tr$status,
      character(1)
    )
  }
  n_total <- length(statuses)
  n_success <- sum(statuses == "success", na.rm = TRUE)
  n_failed <- sum(statuses == "failed", na.rm = TRUE)
  n_remaining <- sum(is.na(statuses) | statuses %in% TASK_RESUMABLE_STATUSES)
  elapsed <- round(result$timing$total, 1)

  if (n_remaining == 0L) {
    if (n_failed == 0L) {
      cli::cli_alert_success(
        "Simulation complete: {n_success}/{n_total} tasks succeeded in {elapsed}s"
      )
    } else {
      cli::cli_alert_warning(
        "Simulation finished: {n_success}/{n_total} tasks succeeded, {n_failed} failed in {elapsed}s"
      )
    }
  } else {
    reasons <- if (
      is.data.frame(grid) &&
        nrow(grid) > 0L &&
        "stop_reason" %in% names(grid)
    ) {
      unique(stats::na.omit(grid$stop_reason))
    } else {
      character()
    }
    reason <- if (length(reasons) > 0L) {
      format_stop_reasons(reasons)
    } else {
      "stopped early"
    }
    cli::cli_alert_warning(
      "Simulation stopped early ({reason}): {n_success} succeeded, {n_failed} failed, {n_remaining} of {n_total} tasks not run"
    )
  }

  if (!is.null(result$checkpoint_path)) {
    path <- normalizePath(
      result$checkpoint_path,
      winslash = "/",
      mustWork = FALSE
    )
    cli::cli_text("Results: {path}")
    if (n_remaining > 0L) {
      # verbatim: the guidance lines are pre-formatted (own indentation, no
      # cli markup), and paths may contain characters cli would interpolate.
      for (line in resume_guidance_lines(path)) {
        cli::cli_verbatim(line)
      }
    }
  }

  print_failure_summary(result)
  invisible(NULL)
}

# F3: as_tibble for simulation results ------------------------------------

# Registered as an S3 method for tibble::as_tibble via registerS3method on load
# (see zzz.R). Returns the per-task summary tibble so tidyverse users can pipe
# results directly: result |> dplyr::filter(...) etc.
as_tibble.bayesim_simulation_result <- function(x, ...) {
  if (!is.null(x$summary) && is.data.frame(x$summary)) {
    tibble::as_tibble(as.data.frame(x$summary))
  } else {
    tibble::tibble()
  }
}
