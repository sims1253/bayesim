# Analysis & reporting layer ---------------------------------------------
#
# Post-simulation summaries, SBC diagnostics, and plotting functions. These
# operate on the `summary` tibble of a `bayesim_simulation_result` (or the
# `task_results` list) and return tidy tibbles or ggplot objects.
#
# MCSE formulas follow rsimsum (Gasparini, 2018): bias MCSE is the sd of the
# estimation errors (est - truth) over sqrt(n), valid under fixed and varying
# truth; coverage MCSE = sqrt(p(1-p)/n); rmse MCSE via the delta method on the
# squared-error mean and variance.

validate_analysis_columns <- function(
  data,
  columns,
  invalid_message,
  unknown_prefix,
  allow_null = TRUE
) {
  if (is.null(columns) && isTRUE(allow_null)) {
    return(invisible(NULL))
  }
  if (!is.character(columns) || anyNA(columns) || !all(nzchar(columns))) {
    stop(bayesim_config_error(invalid_message))
  }
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    stop(bayesim_config_error(paste0(
      unknown_prefix,
      paste(missing, collapse = ", ")
    )))
  }
  invisible(columns)
}

# summarize_simulation ----------------------------------------------------

#' Aggregate simulation results per condition
#'
#' @description Computes per-condition aggregates of the wide summary tibble:
#'   mean and median of each metric column, Monte Carlo standard errors
#'   (MCSE), replicate counts, and failure/convergence-failure rates. Returns a
#'   tidy tibble with one row per condition.
#'
#'   Aggregation follows each metric's declared `summary_type` (E4; see
#'   [Metric]): `"mean"` columns get a `sd / sqrt(n)` MCSE, `"proportion"`
#'   columns (e.g. coverage) get `sqrt(p(1-p) / n)`, and `"none"` columns
#'   (e.g. SBC ranks) are excluded from aggregation. Columns from unknown or
#'   user-defined sources default to `"mean"`. MCSE formulas follow rsimsum
#'   (Gasparini, 2018).
#'
#'   Wide summaries: several metrics legitimately flatten to dozens of columns
#'   each, so the default aggregation can return 100+ columns. Nothing is ever
#'   dropped or truncated. Narrow the output with the `metrics` argument, and
#'   discover a single metric's flattened columns with [metric_cols()]. In
#'   interactive sessions only, a wide default call prints a one-line hint
#'   pointing at these; programmatic and noninteractive use is always silent.
#'
#' @param result A `bayesim_simulation_result` object (uses `result$summary`),
#'   or a data.frame/tibble of per-task metrics. Passing the full result is
#'   preferred: it carries the metrics' `summary_type` declarations.
#' @param by Character vector of grouping columns (conditions). Defaults to the
#'   `data_*`/`fit_*` grid columns found in the summary plus any other
#'   non-numeric condition columns (excluding `task_id` and `status`).
#' @param metrics Character vector of metric columns to aggregate. Defaults to
#'   all numeric columns not in `by` and not metadata
#'   (`task_id`, `rep_idx`, `status`, `*timing*`).
#'
#' @return A tibble with one row per condition: the `by` columns, then for each
#'   metric `<m>_n_used`, `<m>_mean`, `<m>_median`, `<m>_sd`, `<m>_mcse`, plus
#'   `n_reps`, `n_failed`, `failure_rate`. `<m>_n_used` is the number of finite
#'   values used for that metric; failed or non-finite metric values do not
#'   contribute to its aggregate.
#' @export
#' @examples
#' \dontrun{
#' result <- run_simulation(config, progress = FALSE)
#' summarize_simulation(result)
#' summarize_simulation(result, by = "model", metrics = "rmse__value")
#' }
summarize_simulation <- function(result, by = NULL, metrics = NULL) {
  df <- if (is.data.frame(result)) result else result$summary
  if (is.null(df) || !is.data.frame(df)) {
    stop(bayesim_config_error(
      "summarize_simulation requires a simulation result or summary data.frame"
    ))
  }

  # E4: per-metric summary_type declared on the Metric objects (recorded in the
  # result by run_simulation). Column prefix before the first "__" is the
  # metric name. Unknown/user columns default to "mean"; a "coverage"-prefixed
  # column defaults to "proportion" so bare data.frame input still gets the
  # right MCSE.
  summary_types <- if (!is.data.frame(result)) {
    result$metric_summary_types %||% list()
  } else {
    list()
  }
  field_schemas <- if (!is.data.frame(result)) {
    result$metric_field_metadata %||% list()
  } else {
    list()
  }
  field_metadata_for <- function(col) {
    metric_name <- sub("__.*$", "", col)
    suffix <- substring(col, nchar(metric_name) + 3L)
    field <- sub("__.*$", "", suffix)
    schema <- field_schemas[[metric_name]]
    if (S7::S7_inherits(schema, Metric)) {
      return(metric_field_metadata(schema, field))
    }
    if (is.list(schema) && field %in% names(schema)) {
      meta <- schema[[field]]
      if (is.character(meta) && length(meta) == 1L) {
        meta <- list(aggregation = meta)
      }
      if (is.list(meta)) {
        if (is.null(meta$aggregation)) {
          meta$aggregation <- "mean"
        }
        if (is.null(meta$mcse)) {
          meta$mcse <- if (identical(meta$aggregation, "proportion")) {
            "binomial"
          } else if (identical(meta$aggregation, "none")) {
            "none"
          } else {
            "sd"
          }
        }
        return(meta)
      }
    }
    list()
  }
  type_for <- function(col) {
    metadata <- field_metadata_for(col)
    if (!is.null(metadata$aggregation)) {
      return(metadata$aggregation)
    }
    metric_name <- sub("__.*$", "", col)
    declared <- summary_types[[metric_name]]
    if (!is.null(declared)) {
      if (is.list(declared)) {
        field <- sub("^[^_]+__", "", col)
        field <- sub("__.*$", "", field)
        declared <- declared[[field]] %||% declared$aggregation
      }
      if (!is.null(declared)) return(declared)
    }
    if (grepl("^coverage(__|$)", col)) "proportion" else "mean"
  }
  role_for <- function(col) {
    declared <- field_metadata_for(col)$role
    if (!is.null(declared)) {
      return(declared)
    }
    # Compatibility fallback for summaries produced before field schemas were
    # persisted. Counts are bookkeeping, never scientific measures.
    if (grepl("(^|__)n_(obs|draws|ranks)$", col)) {
      return("count")
    }
    ""
  }
  mcse_for <- function(col) {
    field_metadata_for(col)$mcse %||%
      if (identical(type_for(col), "proportion")) "binomial" else "sd"
  }

  if (is.null(by)) {
    # Default grouping: the data_grid/fit_grid condition columns (whatever
    # their type), plus any other non-numeric non-metadata columns.
    grid_cols <- grep("^(data_|fit_)", names(df), value = TRUE)
    other <- names(df)[!vapply(df, is.numeric, logical(1))]
    # Flattened metric payloads use `metric__field`; character payloads such as
    # error messages and artifact paths are outputs, never design conditions.
    other <- other[!grepl("__", other, fixed = TRUE)]
    error_cols <- names(df)[grepl(
      "(^error_(class|message)$|__error_(class|message)$)",
      names(df)
    )]
    by <- setdiff(
      unique(c(grid_cols, other)),
      c("task_id", "status", error_cols)
    )
    # Only scalar atomic columns can group (drop list-columns like fit_formula).
    by <- by[vapply(df[, by, drop = FALSE], is.atomic, logical(1))]
  }
  validate_analysis_columns(
    df,
    by,
    "by must be NULL or a character vector of columns",
    "Unknown grouping column(s): "
  )

  meta_cols <- c(
    "task_id",
    "rep_idx",
    "status",
    "timing_total",
    "timing_warmup",
    "timing_sample"
  )
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  # exclude timing/pure-metadata numeric columns from default metrics
  exclude_meta <- numeric_cols[grepl("^(timing|n_)|(_timing)$", numeric_cols)]
  metric_cols <- setdiff(numeric_cols, c(meta_cols, exclude_meta, by))
  # summary_type = "none" opts a metric out of aggregation (e.g. SBC ranks).
  metric_cols <- metric_cols[vapply(
    metric_cols,
    function(m) type_for(m) != "none" && role_for(m) != "count",
    logical(1)
  )]

  if (!is.null(metrics)) {
    validate_analysis_columns(
      df,
      metrics,
      "metrics must be a character vector of columns",
      "Unknown metric column(s): ",
      allow_null = FALSE
    )
    metric_cols <- intersect(metric_cols, metrics)
  }
  if (length(metric_cols) == 0L) {
    warning("No metric columns found to summarize")
    return(tibble::as_tibble(df[0, , drop = FALSE]))
  }
  # Wide-summary discoverability (interactive only): point at the existing
  # metrics/metric_cols() narrowing instead of ever dropping columns.
  maybe_wide_summary_hint(length(metric_cols), metrics)

  # Ensure status exists for failure-rate computation.
  has_status <- "status" %in% names(df)

  # Split-apply-combine via dplyr if available, else base R.
  groups <- group_ids(df, by)

  rows <- lapply(split(seq_len(nrow(df)), groups, drop = TRUE), function(idx) {
    sub <- df[idx, , drop = FALSE]
    out <- if (length(by) > 0L) as.list(sub[1, by, drop = FALSE]) else list()
    out$n_reps <- length(idx)
    if (has_status) {
      n_fail <- sum(is.na(sub$status) | as.character(sub$status) != "success")
      out$n_failed <- n_fail
      out$failure_rate <- n_fail / length(idx)
    } else {
      out$n_failed <- NA_integer_
      out$failure_rate <- NA_real_
    }
    for (m in metric_cols) {
      vals <- sub[[m]]
      if (has_status) {
        vals <- vals[!is.na(sub$status) & as.character(sub$status) == "success"]
      }
      vals <- vals[is.finite(vals)]
      n <- length(vals)
      out[[paste0(m, "_n_used")]] <- n
      out[[paste0(m, "_mean")]] <- if (n > 0L) mean(vals) else NA_real_
      out[[paste0(m, "_median")]] <- if (n > 0L) {
        stats::median(vals)
      } else {
        NA_real_
      }
      sd_v <- if (n > 1L) stats::sd(vals) else NA_real_
      out[[paste0(m, "_sd")]] <- sd_v
      # E4: MCSE by declared summary_type. "proportion" (coverage-style)
      # columns use sqrt(p(1-p)/n); everything else (including pred-RMSE
      # means) uses the plain sd/sqrt(n) MCSE for what is reported (the mean
      # of per-task values).
      if (n > 1L) {
        if (identical(mcse_for(m), "binomial")) {
          p <- mean(vals)
          out[[paste0(m, "_mcse")]] <- if (is.finite(p) && p >= 0 && p <= 1) {
            sqrt(p * (1 - p) / n)
          } else {
            NA_real_
          }
        } else if (identical(mcse_for(m), "none")) {
          out[[paste0(m, "_mcse")]] <- NA_real_
        } else {
          out[[paste0(m, "_mcse")]] <- sd_v / sqrt(n)
        }
      } else {
        out[[paste0(m, "_mcse")]] <- NA_real_
      }
    }
    out
  })

  # Build a 1-row data.frame per group, then rbind. Uniform schema: rows from
  # different groups always share keys, but we defensively pad any missing
  # column (logical NA, which coerces upward safely in rbind) and align the
  # column order so mixed-type columns never reorder or coerce between groups.
  row_dfs <- lapply(
    rows,
    as.data.frame,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  row_names <- unique(unlist(lapply(row_dfs, names), use.names = FALSE))
  row_dfs <- lapply(row_dfs, function(row) {
    missing <- setdiff(row_names, names(row))
    if (length(missing) > 0L) {
      row[missing] <- rep(list(NA), length(missing))
    }
    row[row_names]
  })

  tibble::as_tibble(do.call(rbind, row_dfs))
}

# Metric-column count at which summarize_simulation() gives a one-line
# interactive hint about narrowing with `metrics = ...` / metric_cols().
# Below this, a wide default summary stays completely silent everywhere.
WIDE_SUMMARY_HINT_THRESHOLD <- 20L

# Interactive-only discoverability hint for wide summaries.
#
# A default summarize_simulation() call over several metrics can legitimately
# aggregate 100+ columns. Rather than dropping or truncating anything, guide
# interactive users toward the existing narrowing capabilities. The hint is
# narrowly triggered: never in noninteractive/programmatic use (tests,
# scripts, reports), never when `metrics` was supplied explicitly, and never
# for narrow summaries.
maybe_wide_summary_hint <- function(
  n_metric_cols,
  metrics = NULL,
  .interactive = interactive()
) {
  if (!isTRUE(.interactive) || !is.null(metrics)) {
    return(invisible(NULL))
  }
  if (
    !is.finite(n_metric_cols) || n_metric_cols <= WIDE_SUMMARY_HINT_THRESHOLD
  ) {
    return(invisible(NULL))
  }
  cli::cli_inform(c(
    i = paste0(
      "Aggregating {n_metric_cols} metric columns. Pass `metrics = <columns>` ",
      "to narrow the output, or use metric_cols() to list one metric's columns ",
      "(see ?summarize_simulation)."
    )
  ))
  invisible(NULL)
}

#' Select flattened metric columns
#'
#' Selects columns belonging to one metric in bayesim's flattened naming
#' scheme, `<metric>__<field>` or `<metric>__<field>__<param>`.
#'
#' @param x A bayesim simulation result (with a `summary` component), or a
#'   per-task summary data frame.
#' @param metric Required character scalar naming the metric prefix.
#' @param fields Optional character vector restricting the field segment.
#' @param as Output form: `"names"` returns a named character vector suitable
#'   for `dplyr::all_of()`; `"long"` returns one row per task and metric value.
#'
#' @return For `as = "names"`, a named character vector whose values are full
#'   column names and whose names are suffixes following the metric prefix. For
#'   `as = "long"`, a tibble with optional `task_id`, `field`, `param`, and
#'   `value` columns.
#' @export
#' @examples
#' fitter <- LinearRegressionFitter()
#' summary <- data.frame(
#'   task_id = "task_1",
#'   posterior_summary__mean__x = 1.2,
#'   posterior_summary__sd__x = 0.3,
#'   check.names = FALSE
#' )
#' metric_cols(summary, "posterior_summary")
#' metric_cols(summary, "posterior_summary", fields = "mean", as = "long")
metric_cols <- function(x, metric, fields = NULL, as = c("names", "long")) {
  as <- match.arg(as)
  if (
    !is.character(metric) ||
      length(metric) != 1L ||
      is.na(metric) ||
      !nzchar(metric)
  ) {
    stop(bayesim_config_error("metric must be a non-empty character scalar"))
  }
  df <- if (is.data.frame(x)) x else x$summary
  if (is.null(df) || !is.data.frame(df)) {
    stop(bayesim_config_error(
      "metric_cols() requires a simulation result or summary data.frame"
    ))
  }

  prefix <- paste0(metric, "__")
  candidates <- names(df)[startsWith(names(df), prefix)]
  suffix <- substring(candidates, nchar(prefix) + 1L)
  parts <- strsplit(suffix, "__", fixed = TRUE)
  field <- vapply(parts, `[[`, character(1), 1L)
  if (!is.null(fields)) {
    keep <- field %in% fields
    candidates <- candidates[keep]
    suffix <- suffix[keep]
    parts <- parts[keep]
    field <- field[keep]
  }
  if (length(candidates) == 0L) {
    flattened <- names(df)[grepl("__", names(df), fixed = TRUE)]
    existing <- sort(unique(sub("__.*$", "", flattened)))
    existing_text <- if (length(existing)) {
      paste(existing, collapse = ", ")
    } else {
      "none"
    }
    stop(bayesim_config_error(paste0(
      "No columns found for metric '",
      metric,
      "'. Metric prefixes present: ",
      existing_text
    )))
  }

  if (as == "names") {
    return(stats::setNames(candidates, suffix))
  }

  rows <- lapply(seq_along(candidates), function(i) {
    param <- if (length(parts[[i]]) > 1L) {
      paste(parts[[i]][-1L], collapse = "__")
    } else {
      NA_character_
    }
    out <- tibble::tibble(
      field = field[[i]],
      param = param,
      value = df[[candidates[[i]]]]
    )
    if ("task_id" %in% names(df)) {
      out <- tibble::add_column(out, task_id = df$task_id, .before = 1L)
    }
    out
  })
  do.call(rbind, rows)
}

# SBC diagnostics ---------------------------------------------------------

#' Extract SBC ranks from a simulation result
#'
#' @description Collects the per-task `rank__by_param` entries (from
#'   [rank_metric()]) into a long tibble with columns `task_id`, `param`,
#'   `rank`, `n_draws`, plus atomic `data_*` and `fit_*` condition columns.
#'   Returns an empty tibble if no rank metric was computed.
#'
#' @param result A `bayesim_simulation_result`.
#' @return A tibble with columns `task_id`, `param`, `rank`, `n_draws`, and
#'   available condition columns.
#' @export
#' @examples
#' \dontrun{
#' result <- run_simulation(config, progress = FALSE)
#' ranks <- sbc_ranks(result)
#' }
sbc_ranks <- function(result) {
  df <- if (is.data.frame(result)) result else result$summary
  empty <- tibble::tibble(
    task_id = character(0),
    param = character(0),
    rank = integer(0),
    n_draws = integer(0),
    n_ranks = integer(0)
  )
  if (is.null(df)) {
    return(empty)
  }
  # rank__by_param__<param> columns hold the per-task per-parameter rank counts.
  # For the single-parameter case the flattener collapses the length-1 named
  # vector to a scalar `rank__by_param` (no __<param> suffix); handle both.
  rank_cols_multi <- grep("^rank__by_param__", names(df), value = TRUE)
  rank_cols_single <- if ("rank__by_param" %in% names(df)) {
    "rank__by_param"
  } else {
    character(0)
  }
  rank_cols <- c(rank_cols_multi, rank_cols_single)
  n_draws_col <- grep("^rank__n_draws$", names(df), value = TRUE)
  if (length(rank_cols) == 0L) {
    return(empty)
  }
  n_draws <- if (length(n_draws_col) == 1L) df[[n_draws_col]] else NA_integer_
  condition_cols <- grep("^(data_|fit_)", names(df), value = TRUE)
  condition_cols <- condition_cols[!grepl("__", condition_cols, fixed = TRUE)]
  condition_cols <- condition_cols[vapply(
    df[, condition_cols, drop = FALSE],
    is.atomic,
    logical(1)
  )]
  rows <- list()
  for (col in rank_cols) {
    param <- sub("^rank__by_param__?", "", col)
    if (param == "") {
      param <- "(single)"
    }
    # Per-variable n_ranks when present (F4); fall back to n_draws + 1.
    n_ranks_col <- paste0("rank__n_ranks__", param)
    n_ranks <- if (n_ranks_col %in% names(df)) {
      as.integer(df[[n_ranks_col]])
    } else {
      NA_integer_
    }
    rank_row <- tibble::tibble(
      task_id = df$task_id,
      param = param,
      rank = as.integer(df[[col]]),
      n_draws = n_draws,
      n_ranks = n_ranks
    )
    if (length(condition_cols) > 0L) {
      rank_row <- cbind(rank_row, df[, condition_cols, drop = FALSE])
    }
    rows[[param]] <- rank_row
  }
  do.call(rbind, rows)
}

#' Plot SBC rank histograms
#'
#' @description One histogram panel per parameter. Requires ggplot2.
#' @param ranks A tibble from [sbc_ranks()], or a `bayesim_simulation_result`.
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' plot_rank_hist(sbc_ranks(result))
#' }
plot_rank_hist <- function(ranks) {
  rlang::check_installed("ggplot2", "to use plot_rank_hist()")
  if (inherits(ranks, "bayesim_simulation_result")) {
    ranks <- sbc_ranks(ranks)
  }
  if (nrow(ranks) == 0L) {
    stop(bayesim_config_error("No rank data; was rank_metric() used?"))
  }
  ggplot2::ggplot(ranks, ggplot2::aes(.data$rank)) +
    ggplot2::geom_histogram(bins = 20, color = "white") +
    ggplot2::facet_wrap(~param, scales = "free_x") +
    ggplot2::labs(x = "rank", y = "count", title = "SBC rank histograms") +
    ggplot2::theme_minimal()
}

#' Plot SBC rank ECDF with simultaneous confidence band
#'
#' @description Plots the empirical CDF of SBC ranks against the uniform CDF
#'   (the diagonal), with a simultaneous confidence band following
#'   Säilynoja, Bürkner, and Vehtari (2022). The band is calibrated so that,
#'   under correct calibration, the *entire* ECDF stays within it with
#'   probability alpha; deviations anywhere along the band therefore indicate
#'   miscalibration at level 1 - alpha.
#'
#'   Ranks are normalized per task: each task's ranks are scaled by that
#'   task's own support, `(rank + 0.5) / n_ranks` with `n_ranks` = support +
#'   1 (kept post-thinning draws + 1). When tasks in a panel have different
#'   supports (e.g. `rank_metric()` with `thin = "auto"` under
#'   autocorrelation), pooling on the panel maximum would squash
#'   small-support tasks' ranks toward zero and manufacture apparent
#'   miscalibration; per-task normalization avoids that artifact, and a
#'   warning notes that the simultaneous band, which assumes iid ranks on a
#'   common support, is then approximate. Legacy results without a recorded
#'   `n_ranks` fall back per task to `n_draws`; with a historical thinning
#'   stride > 1 the true support is unknown, so that fallback only bounds it.
#' @param ranks A tibble from [sbc_ranks()], or a `bayesim_simulation_result`.
#' @param alpha Coverage level of the simultaneous confidence band
#'   (default 0.95).
#' @param by Optional character vector of condition columns to facet by. These
#'   columns are preserved by [sbc_ranks()] for simulation results. Using `by`
#'   computes a separate ECDF and simultaneous band per condition cell instead
#'   of pooling ranks across cells.
#' @references Säilynoja T, Bürkner PC, Vehtari A (2022). Graphical test for
#'   discrete uniformity and its applications in goodness-of-fit evaluation.
#'   *Statistics and Computing*, 32(2).
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' plot_rank_ecdf(sbc_ranks(result))
#' plot_rank_ecdf(sbc_ranks(result), alpha = 0.99)
#' plot_rank_ecdf(sbc_ranks(result), by = "data_n")
#' }
plot_rank_ecdf <- function(ranks, alpha = 0.95, by = NULL) {
  rlang::check_installed("ggplot2", "to use plot_rank_ecdf()")
  if (inherits(ranks, "bayesim_simulation_result")) {
    ranks <- sbc_ranks(ranks)
  }
  if (nrow(ranks) == 0L) {
    stop(bayesim_config_error("No rank data; was rank_metric() used?"))
  }
  if (
    !is.numeric(alpha) ||
      length(alpha) != 1L ||
      is.na(alpha) ||
      alpha <= 0 ||
      alpha >= 1
  ) {
    stop(bayesim_config_error(
      "alpha must be a scalar in (0, 1), got " %+% alpha
    ))
  }
  validate_analysis_columns(
    ranks,
    by,
    "by must be NULL or a character vector",
    "Unknown SBC facet column(s): "
  )

  facet_cols <- c("param", by)
  group_keys <- group_ids(ranks, facet_cols)
  rank_groups <- split(seq_len(nrow(ranks)), group_keys, drop = TRUE)

  heterogeneous_support <- FALSE
  ecdf_data <- do.call(
    rbind,
    lapply(seq_along(rank_groups), function(group_idx) {
      sub <- ranks[rank_groups[[group_idx]], , drop = FALSE]
      n <- nrow(sub)
      # Per-task support S_i (number of possible ranks minus one). Prefer the
      # post-thinning n_ranks (F4); rows without a usable n_ranks fall back
      # to their own n_draws (legacy results: with a historical thinning
      # stride > 1 the true support is unknown and n_draws only bounds it).
      # The panel max rank is the last resort when even n_draws is missing.
      support <- rep(NA_real_, n)
      nr <- if ("n_ranks" %in% names(sub)) as.numeric(sub$n_ranks) else NA_real_
      has_nr <- is.finite(nr) & nr >= 2
      support[has_nr] <- nr[has_nr] - 1
      nd <- if ("n_draws" %in% names(sub)) as.numeric(sub$n_draws) else NA_real_
      has_nd <- !has_nr & is.finite(nd) & nd >= 1
      support[has_nd] <- nd[has_nd]
      if (anyNA(support)) {
        support[is.na(support)] <- max(sub$rank, na.rm = TRUE)
      }
      if (length(unique(nr[is.finite(nr)])) > 1L) {
        heterogeneous_support <<- TRUE
      }
      # Normalized rank on (0,1): rank in 0..S_i maps to (rank+0.5)/(S_i+1),
      # each row by its own support. Sorting the *normalized* values (not the
      # raw ranks) keeps small-support tasks from being squashed toward 0
      # when a panel mixes supports (e.g. thin = "auto" under
      # autocorrelation).
      r <- sort((sub$rank + 0.5) / (support + 1))
      # The band assumes one support per panel; use the panel's largest.
      S <- max(support)
      out <- tibble::tibble(
        .sbc_group = group_idx,
        rank_norm = r,
        ecdf = seq_len(n) / n,
        n = n,
        S = S
      )
      cbind(out, sub[rep(1L, n), facet_cols, drop = FALSE])
    })
  )
  if (heterogeneous_support) {
    cli::cli_warn(c(
      paste0(
        "Some SBC panels mix tasks with different rank supports (n_ranks); ",
        "ranks were normalized per task."
      ),
      "i" = paste0(
        "With heterogeneous support the simultaneous confidence band is ",
        "approximate: it assumes iid ranks uniform on a shared support."
      )
    ))
  }

  # Simultaneous confidence band (Säilynoja et al. 2022) for the ECDF of a
  # uniform sample of size n, evaluated at K = min(n, S + 1) points.
  band_data <- do.call(
    rbind,
    lapply(seq_along(rank_groups), function(group_idx) {
      sub <- ecdf_data[
        ecdf_data$.sbc_group == group_idx,
        ,
        drop = FALSE
      ]
      n <- sub$n[1L]
      K <- max(1L, min(n, sub$S[1L] + 1L))
      band <- sbc_band(n, K = K, conf_level = alpha)
      out <- tibble::tibble(
        .sbc_group = group_idx,
        x = band$x,
        lower = pmax(0, band$lower),
        upper = pmin(1, band$upper)
      )
      cbind(
        out,
        sub[rep(1L, nrow(out)), facet_cols, drop = FALSE]
      )
    })
  )

  ggplot2::ggplot(ecdf_data, ggplot2::aes(.data$rank_norm)) +
    ggplot2::geom_ribbon(
      data = band_data,
      ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper),
      alpha = 0.2,
      fill = "grey70",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "red"
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$ecdf)) +
    ggplot2::facet_wrap(stats::as.formula(
      paste("~", paste(facet_cols, collapse = " + "))
    )) +
    ggplot2::labs(
      x = "normalized rank",
      y = "ECDF",
      title = paste0(
        "SBC rank ECDF with ",
        round(alpha * 100),
        "% simultaneous confidence band"
      )
    ) +
    ggplot2::theme_minimal()
}

# Resolve the `estimand` argument against its legacy `var` alias.
#
# Terminology consistency: performance_measures() uses `estimand` while
# plot_recovery() historically used `var`. `estimand` is now the primary
# name on both; `var` remains a silent compatibility alias. Behavior:
# - neither provided -> error naming the `estimand` argument;
# - both provided and equal -> the shared value (no warning; the alias is
#   compatibility, not deprecation noise);
# - both provided and conflicting -> error (silently picking one would hide
#   a probable copy-paste mistake).
resolve_estimand_alias <- function(estimand, var, fn) {
  if (is.null(estimand) && is.null(var)) {
    stop(bayesim_config_error(paste0(
      fn,
      "() requires a parameter name via `estimand` ",
      "(or the legacy `var` alias)"
    )))
  }
  if (is.null(estimand)) {
    estimand <- var
  }
  if (!is.null(var) && !identical(var, estimand)) {
    stop(bayesim_config_error(paste0(
      fn,
      "(): `estimand` (",
      estimand,
      ") and the legacy alias `var` (",
      var,
      ") disagree; supply only `estimand`"
    )))
  }
  if (
    !is.character(estimand) ||
      length(estimand) != 1L ||
      is.na(estimand) ||
      !nzchar(estimand)
  ) {
    stop(bayesim_config_error(paste0(
      fn,
      "(): estimand must be a single non-empty character string"
    )))
  }
  estimand
}

#' Plot parameter recovery (truth vs posterior estimate)
#'
#' @description Scatter of posterior-mean estimates against true parameter
#'   values, per task, with credible-interval segments. Faceted by a condition
#'   column when `by` is supplied. Requires `posterior_summary_metric` to have
#'   been computed and the truth recorded (E1).
#' @param result A `bayesim_simulation_result`.
#' @param estimand Parameter name (a `vars_of_interest` entry); the preferred
#'   terminology, matching [performance_measures()].
#' @param by Optional name of a condition column to facet by (E7).
#' @param var Legacy alias for `estimand`, kept for backward compatibility
#'   with earlier bayesim versions. It is a silent compatibility alias, not a
#'   deprecated argument, and emits no deprecation warning. When both
#'   `estimand` and `var` are supplied they must name the same parameter;
#'   conflicting values are an error.
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' plot_recovery(result, estimand = "b_x")
#' plot_recovery(result, estimand = "b_x", by = "data_n")
#' # Legacy spelling, identical result:
#' plot_recovery(result, var = "b_x")
#' }
plot_recovery <- function(result, estimand = NULL, by = NULL, var = NULL) {
  rlang::check_installed("ggplot2", "to use plot_recovery()")
  estimand <- resolve_estimand_alias(estimand, var, "plot_recovery")
  df <- result$summary
  if (!is.null(by)) {
    if (
      !is.character(by) ||
        anyNA(by) ||
        !all(nzchar(by)) ||
        !all(by %in% names(df))
    ) {
      stop(bayesim_config_error(
        "by must name existing columns in the simulation summary"
      ))
    }
  }
  mean_col <- paste0("posterior_summary__mean__", estimand)
  lower_col <- paste0("posterior_summary__q_lower__", estimand)
  upper_col <- paste0("posterior_summary__q_upper__", estimand)
  if (!(mean_col %in% names(df))) {
    stop(bayesim_config_error(
      "Posterior summary column " %+%
        mean_col %+%
        " not found; " %+%
        "compute posterior_summary_metric() first."
    ))
  }
  # Truth column (E1): prefer truth__<var>, then a legacy true_params__<var>,
  # then a bare `truth` column.
  truth_candidates <- c(
    paste0("truth__", estimand),
    paste0("true_params__", estimand),
    "truth"
  )
  truth_col <- truth_candidates[truth_candidates %in% names(df)][1]
  plot_df <- data.frame(
    truth = if (!is.na(truth_col) && truth_col %in% names(df)) {
      df[[truth_col]]
    } else {
      NA_real_
    },
    estimate = df[[mean_col]],
    lower = if (lower_col %in% names(df)) df[[lower_col]] else NA_real_,
    upper = if (upper_col %in% names(df)) df[[upper_col]] else NA_real_,
    stringsAsFactors = FALSE
  )
  # Attach the facet column when requested (E7).
  if (!is.null(by) && by %in% names(df)) {
    plot_df$.facet <- df[[by]]
  }
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(.data$truth, .data$estimate)) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      color = "grey50",
      linetype = "dashed"
    ) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::labs(
      x = paste0("true ", estimand),
      y = paste0("posterior mean ", estimand),
      title = paste0("Parameter recovery: ", estimand)
    ) +
    ggplot2::theme_minimal()
  # Posterior-interval segments by default (E7).
  if (lower_col %in% names(df)) {
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        alpha = 0.2
      )
  }
  if (!is.null(by) && by %in% names(df)) {
    p <- p + ggplot2::facet_wrap(~ .data$.facet)
  }
  p
}

#' Plot coverage rates per condition/parameter
#'
#' @description Point-range plot of credible-interval coverage with MCSE error
#'   bars (E7 redesign: was a bar plot of a continuous coverage_mean). Coverage
#'   and its MCSE come from [performance_measures()]; each point is a
#'   condition x estimand cell, with a dashed reference line at the nominal
#'   rate. Requires `posterior_summary_metric()` (and recorded truths, E1).
#' @param result A `bayesim_simulation_result`.
#' @param nominal Nominal coverage rate. Defaults to the interval probability
#'   recorded by the metric schema, or 0.95 for legacy results.
#' @param by Character vector of condition columns (passed to
#'   [performance_measures()]).
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' plot_coverage(result)
#' }
plot_coverage <- function(result, nominal = NULL, by = NULL) {
  rlang::check_installed("ggplot2", "to use plot_coverage()")
  if (is.null(nominal)) {
    schemas <- result$metric_field_metadata %||% list()
    df <- if (is.data.frame(result)) result else result$summary
    lower_params <- sub(
      "^posterior_summary__q_lower__",
      "",
      grep("^posterior_summary__q_lower__", names(df), value = TRUE)
    )
    upper_params <- sub(
      "^posterior_summary__q_upper__",
      "",
      grep("^posterior_summary__q_upper__", names(df), value = TRUE)
    )
    uses_posterior_intervals <- length(intersect(lower_params, upper_params)) >
      0L
    metric_name <- if (uses_posterior_intervals) {
      "posterior_summary"
    } else {
      "coverage"
    }
    fields <- if (uses_posterior_intervals) {
      c("q_lower", "q_upper")
    } else {
      "by_param"
    }
    schema <- schemas[[metric_name]]
    if (S7::S7_inherits(schema, Metric)) {
      schema <- tryCatch(schema@schema, error = function(e) list())
    }
    candidates <- vapply(
      fields,
      function(field) {
        candidate <- if (is.list(schema) && is.list(schema[[field]])) {
          schema[[field]]$nominal
        } else {
          NULL
        }
        if (is.null(candidate)) NA_real_ else as.numeric(candidate)
      },
      numeric(1)
    )
    candidates <- unique(candidates[is.finite(candidates)])
    if (length(candidates) > 1L) {
      stop(bayesim_config_error(paste0(
        "Conflicting nominal levels in ",
        metric_name,
        " interval metadata"
      )))
    }
    nominal <- if (length(candidates) == 1L) candidates[[1L]] else 0.95
  }
  if (
    !is.numeric(nominal) ||
      length(nominal) != 1L ||
      is.na(nominal) ||
      !is.finite(nominal) ||
      nominal <= 0 ||
      nominal >= 1
  ) {
    stop(bayesim_config_error(
      "nominal must be a single finite value in (0, 1)"
    ))
  }
  pm <- performance_measures(result, by = by)
  cov <- pm[pm$measure == "coverage", , drop = FALSE]
  if (nrow(cov) == 0L) {
    stop(bayesim_config_error(paste(
      "No coverage rows: performance_measures() needs truth__ and",
      "posterior_summary__ columns (compute posterior_summary_metric(), E1 truths)."
    )))
  }
  condition_cols <- if (!is.null(by)) {
    by
  } else {
    grep("^(data_|fit_)", names(cov), value = TRUE)
  }
  p <- ggplot2::ggplot(cov, ggplot2::aes(.data$estimand, .data$value)) +
    ggplot2::geom_hline(
      yintercept = nominal,
      color = "grey50",
      linetype = "dashed"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data$value - .data$mcse,
        ymax = .data$value + .data$mcse
      ),
      width = 0.2
    ) +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "estimand",
      y = paste0("coverage (nominal ", nominal, ")"),
      title = "Credible-interval coverage with MCSE"
    ) +
    ggplot2::theme_minimal()
  if (length(condition_cols) > 0L) {
    p <- p +
      ggplot2::facet_wrap(stats::as.formula(
        paste("~", paste(condition_cols, collapse = " + "))
      ))
  }
  p
}

#' Plot a metric across conditions
#'
#' @description General-purpose plot of a metric column (or its aggregated mean)
#'   against a condition variable, with optional facets. Requires ggplot2.
#' @param result A `bayesim_simulation_result` or summary tibble.
#' @param metric Name of the metric column to plot.
#' @param x Optional x-axis variable (condition). Defaults to `task_id`.
#' @param facets Optional facet variable.
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' plot_metric(result, "rmse__value", x = "n")
#' }
plot_metric <- function(result, metric, x = NULL, facets = NULL) {
  rlang::check_installed("ggplot2", "to use plot_metric()")
  df <- if (is.data.frame(result)) result else result$summary
  if (!(metric %in% names(df))) {
    stop(bayesim_config_error(
      "Metric column " %+% metric %+% " not found in summary"
    ))
  }
  df$.metric <- df[[metric]]
  if (is.null(x)) {
    x <- "task_id"
  }
  if (
    !is.character(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !nzchar(x) ||
      !x %in% names(df)
  ) {
    stop(bayesim_config_error(
      "x must name an existing column in the summary"
    ))
  }
  if (!is.null(facets)) {
    if (
      !is.character(facets) ||
        anyNA(facets) ||
        !all(nzchar(facets)) ||
        !all(facets %in% names(df))
    ) {
      stop(bayesim_config_error(
        "facets must name existing columns in the summary"
      ))
    }
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(.data[[x]], .data$.metric)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::labs(y = metric, title = metric) +
    ggplot2::theme_minimal()
  if (!is.null(facets)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facets)))
  }
  p
}

# performance_measures (Morris et al. 2019 layer) -------------------------

#' Simulation-method performance measures with Monte-Carlo standard errors
#'
#' @description
#' For condition cells with fixed data-generating truth, computes the Morris,
#' White & Crowther (2019, *Stat Med*) estimator-performance measures:
#' **bias**, **empirical SE**, **MSE**, **coverage**, **average model SE**, and
#' `n_sim`, each with its MCSE. When truth varies across replicates (as in a
#' prior-predictive SBC study), the fixed-truth Morris names are not returned.
#' Instead, the error distribution is described by `mean_error`, `error_sd`,
#' and `error_mse`; coverage and average model SE remain available.
#'
#' For each parameter the function pairs the data-generating `truth__<param>`
#' column with the per-task `posterior_summary__*__<param>` columns (point
#' estimate `mean`, posterior `sd`, and interval `q_lower`/`q_upper`). Coverage
#' uses the interval when available; otherwise it falls back to a
#' `coverage__by_param__<param>` column if present.
#'
#' Fixed-truth MCSE formulas follow Morris et al. / rsimsum:
#' \itemize{
#'   \item bias MCSE = sd(est - truth) / sqrt(n)
#'   \item empSE MCSE = sd / sqrt(2(n-1))
#'   \item MSE MCSE = sqrt(Var((est-truth)^2) / n)
#'   \item coverage MCSE = sqrt(p(1-p) / n)
#'   \item modelSE MCSE = sd(posterior_sd) / sqrt(n)
#' }
#' The bias MCSE uses the sd of the estimation errors `est - truth`, which is
#' valid under fixed and varying truth alike (with fixed truth the truth is
#' constant, so sd(est - truth) = sd(est)). For varying truth, `mean_error`
#' uses `sd(est-truth) / sqrt(n)`, `error_sd` uses `error_sd / sqrt(2(n-1))`,
#' and `error_mse` uses `sqrt(Var((est-truth)^2) / n)`.
#'
#' @param result A `bayesim_simulation_result` (uses `$summary`), or a
#'   data.frame of per-task metrics.
#' @param estimand Optional character; a single parameter name. When NULL
#'   (default), all parameters with both a `truth__*` and
#'   `posterior_summary__mean__*` column are analyzed.
#' @param estimator Character scalar naming the per-task point estimate to use:
#'   one of `"mean"` (default), `"median"`. Selects the corresponding
#'   `posterior_summary__<estimator>__<param>` column.
#' @param by Character vector of condition/grouping columns. Defaults to the
#'   data_grid and fit_grid columns found in the summary (excluding task_id,
#'   rep_idx, status, and metric columns).
#'
#' @return A tidy tibble with columns: the `by` columns, `estimand`,
#'   `measure`, `value`, `mcse`, `n_sim`, `truth_mode`. One row per condition x
#'   estimand x measure. Fixed-truth measures use `bias`, `emp_se`, and `mse`;
#'   varying-truth measures use `mean_error`, `error_sd`, and `error_mse`.
#'   Both modes may also include `coverage`, `model_se`, and `n_sim`.
#' @export
#' @references Morris, White & Crowther (2019), *Using simulation studies to
#'   evaluate statistical methods*, Statistics in Medicine.
#' @examples
#' \dontrun{
#' result <- run_simulation(config, progress = FALSE)
#' performance_measures(result)
#' performance_measures(result, estimand = "x", by = "model")
#' }
performance_measures <- function(
  result,
  estimand = NULL,
  estimator = c("mean", "median"),
  by = NULL
) {
  estimator <- match.arg(estimator)
  df <- if (is.data.frame(result)) result else result$summary
  if (is.null(df) || !is.data.frame(df)) {
    stop(bayesim_config_error(
      "performance_measures requires a simulation result or summary data.frame"
    ))
  }

  # Default grouping: the data_grid/fit_grid condition columns the summary was
  # enriched with (data_*/fit_* prefixes, including fit_model for
  # model-comparison studies). Metric columns carry "__" and are excluded.
  if (is.null(by)) {
    by <- grep("^(data_|fit_)", names(df), value = TRUE)
    by <- by[!grepl("__", by, fixed = TRUE)]
    # Keep only atomic condition columns (drop e.g. fit_formula list-columns).
    by <- by[vapply(
      df[, by, drop = FALSE],
      function(col) {
        is.numeric(col) ||
          is.character(col) ||
          is.factor(col) ||
          is.logical(col)
      },
      logical(1)
    )]
  }
  validate_analysis_columns(
    df,
    by,
    "by must be NULL or a character vector of columns",
    "Unknown grouping column(s): ",
    allow_null = FALSE
  )

  # Discover estimands: parameters with both truth__ and posterior_summary__mean__.
  truth_cols <- grep("^truth__", names(df), value = TRUE)
  est_col_tmpl <- paste0("posterior_summary__", estimator, "__")
  available_estimands <- vapply(
    truth_cols,
    function(tc) {
      param <- sub("^truth__", "", tc)
      paste0(est_col_tmpl, param) %in% names(df)
    },
    logical(1)
  )
  estimands <- if (is.null(estimand)) {
    sub("^truth__", "", truth_cols[available_estimands])
  } else {
    estimand
  }
  if (length(estimands) == 0L) {
    stop(bayesim_config_error(paste(
      "No estimands found: performance_measures needs both truth__<param> and",
      paste0("posterior_summary__", estimator, "__<param>"),
      "columns. Compute posterior_summary_metric() and record truths (E1)."
    )))
  }

  rows <- list()
  for (param in estimands) {
    truth_col <- paste0("truth__", param)
    est_col <- paste0(est_col_tmpl, param)
    sd_col <- paste0("posterior_summary__sd__", param)
    lo_col <- paste0("posterior_summary__q_lower__", param)
    hi_col <- paste0("posterior_summary__q_upper__", param)
    cov_col <- paste0("coverage__by_param__", param)

    if (!truth_col %in% names(df) || !est_col %in% names(df)) {
      next
    }

    split_cols <- if (length(by)) {
      df[, by, drop = FALSE]
    } else {
      data.frame(.group = rep(1, nrow(df)))
    }
    splits <- group_ids(split_cols, names(split_cols))

    for (lev in unique(splits)) {
      sel <- which(splits == lev)
      if ("status" %in% names(df)) {
        sel <- sel[
          !is.na(df$status[sel]) & as.character(df$status[sel]) == "success"
        ]
      }
      if (length(sel) == 0L) {
        next
      }
      est <- df[[est_col]][sel]
      truth <- df[[truth_col]][sel]
      ok <- is.finite(est) & is.finite(truth)
      est <- est[ok]
      truth <- truth[ok]
      n <- length(est)
      if (n == 0L) {
        next
      }
      errs <- est - truth

      # Condition cell values for the `by` columns (empty when no conditions).
      cond <- if (length(by)) as.list(df[sel[1L], by, drop = FALSE]) else list()
      truth_varies <- length(unique(truth)) > 1L
      truth_mode <- if (truth_varies) "varying" else "fixed"

      add <- function(measure, value, mcse) {
        row <- c(
          cond,
          list(
            estimand = param,
            measure = measure,
            value = value,
            mcse = mcse,
            n_sim = n,
            truth_mode = truth_mode
          )
        )
        rows[[length(rows) + 1L]] <<- row
      }

      # The bias MCSE uses the sd of the estimation errors in both truth
      # modes: with fixed truth sd(est - truth) = sd(est), and when truth
      # varies (prior-predictive or otherwise) the replicate-level error is
      # the only valid reference for bias uncertainty — using sd(est) would
      # mostly measure variation in the data-generating truth rather than
      # simulation error. The empirical-spread measures follow the same
      # split (emp_se from est under fixed truth, error_sd from errs
      # otherwise).
      error_sd <- if (n > 1L) stats::sd(errs) else NA_real_
      sq <- errs^2
      if (truth_varies) {
        add(
          "mean_error",
          mean(errs),
          if (n > 1L) error_sd / sqrt(n) else NA_real_
        )
        add(
          "error_sd",
          error_sd,
          if (n > 2L) error_sd / sqrt(2 * (n - 1L)) else NA_real_
        )
        add(
          "error_mse",
          mean(sq),
          if (n > 1L) sqrt(stats::var(sq) / n) else NA_real_
        )
      } else {
        add(
          "bias",
          mean(errs),
          if (n > 1L) stats::sd(errs) / sqrt(n) else NA_real_
        )
        emp_se <- if (n > 1L) stats::sd(est) else NA_real_
        add(
          "emp_se",
          emp_se,
          if (n > 2L) emp_se / sqrt(2 * (n - 1L)) else NA_real_
        )
        add(
          "mse",
          mean(sq),
          if (n > 1L) sqrt(stats::var(sq) / n) else NA_real_
        )
      }
      # average model SE = mean(posterior_sd); MCSE = sd(posterior_sd)/sqrt(n)
      if (sd_col %in% names(df)) {
        pse <- df[[sd_col]][sel][ok]
        pse <- pse[is.finite(pse)]
        add(
          "model_se",
          if (length(pse) > 0L) mean(pse) else NA_real_,
          if (length(pse) > 1L) {
            stats::sd(pse) / sqrt(length(pse))
          } else {
            NA_real_
          }
        )
      }
      # coverage = P(truth in [lo, hi]); MCSE = sqrt(p(1-p)/n)
      if (lo_col %in% names(df) && hi_col %in% names(df)) {
        lo <- df[[lo_col]][sel][ok]
        hi <- df[[hi_col]][sel][ok]
        covered <- (truth >= lo) & (truth <= hi)
        covered <- covered[!is.na(covered)]
        p <- if (length(covered)) mean(covered) else NA_real_
        ncv <- length(covered)
        add("coverage", p, if (ncv > 1L) sqrt(p * (1 - p) / ncv) else NA_real_)
      } else if (cov_col %in% names(df)) {
        fallback <- df[[cov_col]][sel]
        fallback <- fallback[is.finite(fallback)]
        p <- if (length(fallback)) mean(fallback) else NA_real_
        ncv <- length(fallback)
        add("coverage", p, if (ncv > 1L) sqrt(p * (1 - p) / ncv) else NA_real_)
      }
      add("n_sim", n, NA_real_)
    }
  }

  tibble::as_tibble(do.call(
    rbind,
    lapply(rows, as.data.frame, stringsAsFactors = FALSE)
  ))
}

# n_replicates_for_target (Workstream I5) ----------------------------------

#' Required number of replicates for a target MCSE
#'
#' @description
#' Inverts the Monte-Carlo standard error (MCSE) formulas used by
#' [performance_measures()] to return the number of replicates `n` required to
#' achieve a target MCSE. Useful for planning a simulation study (a precision /
#' power calculation) before running it.
#'
#' The MCSE formulas (Morris et al. 2019; rsimsum) invert to:
#' \itemize{
#'   \item **coverage** (binary 0/1 outcome):
#'     `MCSE = sqrt(p (1 - p) / n)`  =>  `n = p (1 - p) / MCSE^2`.
#'     The default `p = 0.5` is the conservative max-variance case (it
#'     maximises `p (1 - p)` and therefore gives the largest required `n`).
#'   \item **continuous** metrics (bias / mean / model SE):
#'     `MCSE = sd / sqrt(n)`  =>  `n = (sd / MCSE)^2`.
#'     Requires an assumed standard deviation `assumed_sd` for the per-replicate
#'     point estimate (e.g. a guess at the empirical SE of the estimator).
#' }
#'
#' The returned value is `ceiling(n)` so it is always a whole number of
#' replicates.
#'
#' @param target_mcse Numeric scalar > 0. The MCSE you want to achieve.
#' @param metric_type Character scalar: `"coverage"` (binary coverage metrics) or
#'   `"continuous"` (bias / mean / model SE metrics).
#' @param p Numeric scalar in `[0, 1]`. Assumed coverage probability. Only used
#'   when `metric_type = "coverage"`. Defaults to `0.5`, the
#'   variance-maximising (most conservative) choice.
#' @param assumed_sd Numeric scalar > 0. Assumed standard deviation of the
#'   per-replicate point estimate. Required when
#'   `metric_type = "continuous"`; ignored otherwise.
#'
#' @return An integer scalar: the number of replicates required (the ceiling of
#'   the inverted-MCSE `n`).
#' @export
#' @seealso [performance_measures()] for the MCSE formulas being inverted.
#' @examples
#' \dontrun{
#' # Coverage: target MCSE 0.03 under the conservative p = 0.5.
#' n_replicates_for_target(0.03, "coverage")
#'
#' # Coverage at an assumed rate of 0.9.
#' n_replicates_for_target(0.03, "coverage", p = 0.9)
#'
#' # Continuous metric (e.g. bias) with assumed sd of the point estimate.
#' n_replicates_for_target(0.05, "continuous", assumed_sd = 0.5)
#' }
n_replicates_for_target <- function(
  target_mcse,
  metric_type = c("coverage", "continuous"),
  p = 0.5,
  assumed_sd = NULL
) {
  metric_type <- match.arg(metric_type)

  if (
    !is.numeric(target_mcse) ||
      length(target_mcse) != 1L ||
      is.na(target_mcse) ||
      target_mcse <= 0
  ) {
    stop(bayesim_config_error(
      "target_mcse must be a single positive numeric value"
    ))
  }

  if (metric_type == "coverage") {
    if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 0 || p > 1) {
      stop(bayesim_config_error(
        "p must be a single numeric value in [0, 1] for metric_type 'coverage'"
      ))
    }
    n <- p * (1 - p) / target_mcse^2
  } else {
    if (
      is.null(assumed_sd) ||
        !is.numeric(assumed_sd) ||
        length(assumed_sd) != 1L ||
        is.na(assumed_sd) ||
        assumed_sd <= 0
    ) {
      stop(bayesim_config_error(
        "assumed_sd must be a single positive numeric value for metric_type 'continuous'"
      ))
    }
    n <- (assumed_sd / target_mcse)^2
  }

  as.integer(ceiling(n))
}

# report (Workstream I4) ---------------------------------------------------

#' Render a simulation-study report
#'
#' @description
#' Renders a standard Quarto HTML report for a `bayesim_simulation_result`,
#' covering the study design table, [performance_measures()], parameter
#' recovery plots per estimand, credible-interval coverage, and SBC rank ECDF
#' panels (when a rank metric was computed). The report template lives at
#' `inst/report/simulation-report.qmd`.
#'
#' Each section is wrapped in `tryCatch`, so a missing metric (e.g. no rank
#' data, no recorded truths) degrades gracefully instead of failing the whole
#' render.
#'
#' Requires the `quarto` R package AND the Quarto CLI. If the CLI is not
#' available, an informative error is thrown pointing to <https://quarto.org>.
#'
#' `render_report()` was previously named `report()`; the old name collided
#' with the generic of the easystats *report* package and now lives on as a
#' deprecated alias (see [report()]).
#'
#' @param result A `bayesim_simulation_result` from [run_simulation()].
#' @param output_file Path to the rendered HTML output file (default
#'   `"bayesim-report.html"`).
#' @param open Logical scalar. When `TRUE` (the default in interactive
#'   sessions) the rendered report is opened in a viewer/browser. Forwarded to
#'   [quarto::quarto_render()] via its `open` handling where supported; on
#'   systems without that argument the file path is still returned.
#'
#' @return The path to the rendered HTML file (invisibly).
#' @export
#' @seealso [performance_measures()], [plot_recovery()], [plot_coverage()],
#'   [plot_rank_ecdf()].
#' @examples
#' \dontrun{
#' result <- run_simulation(config, progress = FALSE)
#' render_report(result, output_file = "my-study.html")
#' }
render_report <- function(
  result,
  output_file = "bayesim-report.html",
  open = interactive()
) {
  if (!inherits(result, "bayesim_simulation_result")) {
    stop(bayesim_config_error(
      "render_report() requires a bayesim_simulation_result object"
    ))
  }

  rlang::check_installed("quarto", "to render simulation reports")
  if (is.null(quarto::quarto_path())) {
    stop(bayesim_config_error(paste(
      "The Quarto CLI was not found. Install it from https://quarto.org",
      "and ensure it is on your PATH, then call render_report() again."
    )))
  }

  template <- system.file(
    "report",
    "simulation-report.qmd",
    package = "bayesim"
  )
  if (!nzchar(template)) {
    stop(bayesim_config_error(
      "Report template 'inst/report/simulation-report.qmd' is not installed."
    ))
  }

  # Persist the result to a temp rds the template will read via params.
  result_path <- tempfile(fileext = ".rds")
  saveRDS(result, result_path)
  on.exit(unlink(result_path), add = TRUE)

  # Quarto rejects paths in output_file and writes next to the input, so copy
  # the template to a scratch dir, render there, then move the html into place.
  render_dir <- tempfile("bayesim-report-")
  dir.create(render_dir)
  on.exit(unlink(render_dir, recursive = TRUE), add = TRUE)
  render_input <- file.path(render_dir, basename(template))
  file.copy(template, render_input)

  quarto::quarto_render(
    input = render_input,
    output_file = basename(output_file),
    execute_params = list(result_path = result_path)
  )

  rendered <- file.path(render_dir, basename(output_file))
  if (!file.exists(rendered)) {
    stop(bayesim_config_error(
      "Quarto did not produce the expected report file"
    ))
  }
  out <- normalizePath(output_file, mustWork = FALSE)
  file.copy(rendered, out, overwrite = TRUE)
  if (isTRUE(open) && file.exists(out)) {
    tryCatch(utils::browseURL(out), error = function(e) NULL)
  }

  invisible(out)
}

# Once-per-session flag for the report() deprecation warning. Unlike
# .warn_once() in metrics-built-in.R (reset at every run_simulation()), a
# deprecation warning should fire at most once per session, not once per run.
.report_deprecated_env <- new.env(parent = emptyenv())

#' Deprecated alias for render_report()
#'
#' @description
#' `report()` was renamed to [render_report()] because the old name collided
#' with the generic of the easystats *report* package. The alias forwards to
#' `render_report()` and emits a deprecation warning once per session; it will
#' be removed in a future release.
#'
#' @param estimands Ignored. The argument was accepted (and documented as
#'   informational only) by `report()` but never used; it is kept in the
#'   signature purely so existing calls do not error.
#' @inheritParams render_report
#' @export
#' @keywords internal
report <- function(
  result,
  output_file = "bayesim-report.html",
  open = interactive(),
  estimands = NULL
) {
  if (is.null(.report_deprecated_env$warned)) {
    .report_deprecated_env$warned <- TRUE
    cli::cli_warn(c(
      "{.fn report} was renamed to {.fn render_report} and is deprecated.",
      i = "The alias keeps working but will be removed in a future release."
    ))
  }
  render_report(result, output_file = output_file, open = open)
}
