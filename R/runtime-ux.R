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
#' @param pilot Logical; if TRUE, time a single pilot task to estimate total
#'   wall-clock time (default FALSE). Experimental.
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
      if (isTRUE(fitter@supports_loo)) "loo"
    )
  } else {
    character()
  }
  unmet <- setdiff(metric_needs, caps)

  daemons_set <- isTRUE(mirai::daemons_set())
  n_compile <- if (
    inherits(fitter, "BrmsFitter") && !is.null(config@fit_grid)
  ) {
    # Distinct model specs (formula/family combinations) → bank compiles each.
    nrow(config@fit_grid)
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
    n_compile = n_compile
  )

  if (condensed) {
    workers_str <- if (daemons_set) {
      paste0("; ", mirai::daemons()$daemons, " daemons")
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
    if (pilot) cli::cli_text("(pilot timing not implemented in this build)")
  }

  invisible(info)
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
