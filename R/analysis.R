# Analysis & reporting layer ---------------------------------------------
#
# Post-simulation summaries, SBC diagnostics, and plotting functions. These
# operate on the `summary` tibble of a `bayesim_simulation_result` (or the
# `task_results` list) and return tidy tibbles or ggplot objects.
#
# MCSE formulas follow rsimsum (Gasparini, 2018): bias MCSE = sd/sqrt(n);
# coverage MCSE = sqrt(p(1-p)/n); rmse MCSE via the delta method on the
# squared-error mean and variance.

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
#'   metric `<m>_mean`, `<m>_median`, `<m>_sd`, `<m>_mcse`, plus `n_reps`,
#'   `n_failed`, `failure_rate`.
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
    result$metric_summary_types
  } else {
    NULL
  }
  type_for <- function(col) {
    metric_name <- sub("__.*$", "", col)
    declared <- summary_types[[metric_name]]
    if (!is.null(declared)) {
      return(declared)
    }
    if (grepl("^coverage(__|$)", col)) "proportion" else "mean"
  }

  if (is.null(by)) {
    # Default grouping: the data_grid/fit_grid condition columns (whatever
    # their type), plus any other non-numeric non-metadata columns.
    grid_cols <- grep("^(data_|fit_)", names(df), value = TRUE)
    other <- names(df)[!vapply(df, is.numeric, logical(1))]
    by <- setdiff(unique(c(grid_cols, other)), c("task_id", "status"))
    # Only scalar atomic columns can group (drop list-columns like fit_formula).
    by <- by[vapply(df[, by, drop = FALSE], is.atomic, logical(1))]
  }
  by <- intersect(by, names(df))

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
    function(m) type_for(m) != "none",
    logical(1)
  )]

  if (!is.null(metrics)) {
    metric_cols <- intersect(metric_cols, metrics)
  }
  if (length(metric_cols) == 0L) {
    warning("No metric columns found to summarize")
    return(tibble::as_tibble(df[0, , drop = FALSE]))
  }

  # Ensure status exists for failure-rate computation.
  has_status <- "status" %in% names(df)

  # Split-apply-combine via dplyr if available, else base R.
  groups <- if (length(by) > 0L) {
    do.call(paste, c(df[, by, drop = FALSE], sep = "|"))
  } else {
    rep("__all__", nrow(df))
  }

  rows <- lapply(split(seq_len(nrow(df)), groups), function(idx) {
    sub <- df[idx, , drop = FALSE]
    out <- if (length(by) > 0L) as.list(sub[1, by, drop = FALSE]) else list()
    out$n_reps <- length(idx)
    if (has_status) {
      n_fail <- sum(!sub$status %in% "success", na.rm = TRUE)
      out$n_failed <- n_fail
      out$failure_rate <- n_fail / length(idx)
    } else {
      out$n_failed <- NA_integer_
      out$failure_rate <- NA_real_
    }
    for (m in metric_cols) {
      vals <- sub[[m]]
      vals <- vals[is.finite(vals)]
      n <- length(vals)
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
        if (identical(type_for(m), "proportion")) {
          p <- mean(vals)
          out[[paste0(m, "_mcse")]] <- sqrt(p * (1 - p) / n)
        } else {
          out[[paste0(m, "_mcse")]] <- sd_v / sqrt(n)
        }
      } else {
        out[[paste0(m, "_mcse")]] <- NA_real_
      }
    }
    out
  })

  tibble::as_tibble(do.call(rbind, lapply(rows, as.data.frame)))
}

# SBC diagnostics ---------------------------------------------------------

#' Extract SBC ranks from a simulation result
#'
#' @description Collects the per-task `rank__by_param` entries (from
#'   [rank_metric()]) into a long tibble with columns `task_id`, `param`,
#'   `rank`, `n_draws`. Returns an empty tibble if no rank metric was computed.
#'
#' @param result A `bayesim_simulation_result`.
#' @return A tibble with columns `task_id`, `param`, `rank`, `n_draws`.
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
    rows[[param]] <- tibble::tibble(
      task_id = df$task_id,
      param = param,
      rank = as.integer(df[[col]]),
      n_draws = n_draws,
      n_ranks = n_ranks
    )
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
#' @param ranks A tibble from [sbc_ranks()], or a `bayesim_simulation_result`.
#' @param alpha Coverage level of the simultaneous confidence band
#'   (default 0.95).
#' @references Säilynoja T, Bürkner PC, Vehtari A (2022). Graphical test for
#'   discrete uniformity and its applications in goodness-of-fit evaluation.
#'   *Statistics and Computing*, 32(2).
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' plot_rank_ecdf(sbc_ranks(result))
#' plot_rank_ecdf(sbc_ranks(result), alpha = 0.99)
#' }
plot_rank_ecdf <- function(ranks, alpha = 0.95) {
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

  params <- unique(ranks$param)

  ecdf_data <- do.call(
    rbind,
    lapply(params, function(p) {
      sub <- ranks[ranks$param == p, , drop = FALSE]
      n <- nrow(sub)
      # Prefer post-thinning n_ranks (F4); fall back to n_draws for old results.
      S <- NA
      if ("n_ranks" %in% names(sub) && any(is.finite(sub$n_ranks))) {
        S <- max(sub$n_ranks - 1L, na.rm = TRUE)
      }
      if (!is.finite(S) || S < 1L) {
        S <- max(sub$n_draws, na.rm = TRUE)
      }
      if (!is.finite(S) || S < 1L) {
        S <- max(sub$rank, na.rm = TRUE)
      }
      # Normalized rank on [0,1]: rank in 0..S, map to (rank+0.5)/(S+1).
      r <- (sort(sub$rank) + 0.5) / (S + 1)
      tibble::tibble(
        param = p,
        rank_norm = r,
        ecdf = seq_len(n) / n,
        n = n,
        S = S
      )
    })
  )

  # Simultaneous confidence band (Säilynoja et al. 2022) for the ECDF of a
  # uniform sample of size n, evaluated at K = min(n, S + 1) points.
  band_data <- do.call(
    rbind,
    lapply(params, function(p) {
      sub <- ecdf_data[ecdf_data$param == p, , drop = FALSE]
      n <- sub$n[1L]
      K <- max(1L, min(n, sub$S[1L] + 1L))
      band <- sbc_band(n, K = K, conf_level = alpha)
      tibble::tibble(
        param = p,
        x = band$x,
        lower = pmax(0, band$lower),
        upper = pmin(1, band$upper)
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
    ggplot2::facet_wrap(~param) +
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

#' Plot parameter recovery (truth vs posterior estimate)
#'
#' @description Scatter of posterior-mean estimates against true parameter
#'   values, per task, with credible-interval segments. Faceted by a condition
#'   column when `by` is supplied. Requires `posterior_summary_metric` to have
#'   been computed and the truth recorded (E1).
#' @param result A `bayesim_simulation_result`.
#' @param var Parameter name (a vars_of_interest entry).
#' @param by Optional name of a condition column to facet by (E7).
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' plot_recovery(result, "b_x")
#' plot_recovery(result, "b_x", by = "data_n")
#' }
plot_recovery <- function(result, var, by = NULL) {
  rlang::check_installed("ggplot2", "to use plot_recovery()")
  df <- result$summary
  mean_col <- paste0("posterior_summary__mean__", var)
  lower_col <- paste0("posterior_summary__q_lower__", var)
  upper_col <- paste0("posterior_summary__q_upper__", var)
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
    paste0("truth__", var),
    paste0("true_params__", var),
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
      x = paste0("true ", var),
      y = paste0("posterior mean ", var),
      title = paste0("Parameter recovery: ", var)
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
#' @param nominal Nominal coverage rate (default 0.95).
#' @param by Character vector of condition columns (passed to
#'   [performance_measures()]).
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' plot_coverage(result)
#' }
plot_coverage <- function(result, nominal = 0.95, by = NULL) {
  rlang::check_installed("ggplot2", "to use plot_coverage()")
  pm <- performance_measures(result, by = by)
  cov <- pm[pm$measure == "coverage", , drop = FALSE]
  if (nrow(cov) == 0L) {
    stop(bayesim_config_error(paste(
      "No coverage rows: performance_measures() needs truth__ and",
      "posterior_summary__ columns (compute posterior_summary_metric(), E1 truths)."
    )))
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
#' Computes the Morris, White & Crowther (2019, *Stat Med*) estimator-
#' performance measures — **bias**, **empirical SE**, **MSE**, **coverage**,
#' **average model SE**, and `n_sim` — for each estimand (parameter) and
#' condition cell, each with its MCSE. This is the function a methods paper
#' actually needs; it is the centerpiece of the analysis layer (E3).
#'
#' For each parameter the function pairs the data-generating `truth__<param>`
#' column with the per-task `posterior_summary__*__<param>` columns (point
#' estimate `mean`, posterior `sd`, and interval `q_lower`/`q_upper`). Coverage
#' uses the interval when available; otherwise it falls back to a
#' `coverage__by_param__<param>` column if present.
#'
#' MCSE formulas follow Morris et al. / rsimsum:
#' \itemize{
#'   \item bias MCSE = sd(point_est) / sqrt(n)
#'   \item empSE MCSE = sd / sqrt(2(n-1))
#'   \item MSE MCSE = sqrt(Var((est-truth)^2) / n)
#'   \item coverage MCSE = sqrt(p(1-p) / n)
#'   \item modelSE MCSE = sd(posterior_sd) / sqrt(n)
#' }
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
#'   `measure`, `value`, `mcse`, `n_sim`. One row per condition x estimand x
#'   measure. `measure` is one of `bias`, `emp_se`, `mse`, `coverage`,
#'   `model_se`, `n_sim`.
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
    splits <- interaction(split_cols, drop = TRUE)

    for (lev in levels(splits)) {
      sel <- which(splits == lev)
      est <- df[[est_col]][sel]
      truth <- df[[truth_col]][sel]
      ok <- !is.na(est) & !is.na(truth)
      est <- est[ok]
      truth <- truth[ok]
      n <- length(est)
      if (n == 0L) {
        next
      }
      errs <- est - truth

      # Condition cell values for the `by` columns (empty when no conditions).
      cond <- if (length(by)) as.list(df[sel[1L], by, drop = FALSE]) else list()

      add <- function(measure, value, mcse) {
        row <- c(
          cond,
          list(
            estimand = param,
            measure = measure,
            value = value,
            mcse = mcse,
            n_sim = n
          )
        )
        rows[[length(rows) + 1L]] <<- row
      }

      # bias = mean(est - truth); MCSE = sd(est)/sqrt(n)
      add(
        "bias",
        mean(errs),
        if (n > 1L) stats::sd(est) / sqrt(n) else NA_real_
      )
      # empirical SE = sd(est); MCSE = sd/sqrt(2(n-1))
      emp_se <- if (n > 1L) stats::sd(est) else NA_real_
      add(
        "emp_se",
        emp_se,
        if (n > 2L) emp_se / sqrt(2 * (n - 1L)) else NA_real_
      )
      # MSE = mean((est-truth)^2); MCSE = sqrt(Var(err^2)/n)
      sq <- errs^2
      add("mse", mean(sq), if (n > 1L) sqrt(stats::var(sq) / n) else NA_real_)
      # average model SE = mean(posterior_sd); MCSE = sd(posterior_sd)/sqrt(n)
      if (sd_col %in% names(df)) {
        pse <- df[[sd_col]][sel][ok]
        add(
          "model_se",
          mean(pse, na.rm = TRUE),
          if (sum(!is.na(pse)) > 1L) {
            stats::sd(pse, na.rm = TRUE) / sqrt(sum(!is.na(pse)))
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
        p <- mean(covered, na.rm = TRUE)
        ncv <- sum(!is.na(covered))
        add("coverage", p, if (ncv > 1L) sqrt(p * (1 - p) / ncv) else NA_real_)
      } else if (cov_col %in% names(df)) {
        p <- mean(df[[cov_col]][sel], na.rm = TRUE)
        ncv <- sum(!is.na(df[[cov_col]][sel]))
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
#' @param result A `bayesim_simulation_result` from [run_simulation()].
#' @param output_file Path to the rendered HTML output file (default
#'   `"bayesim-report.html"`).
#' @param open Logical scalar. When `TRUE` (the default in interactive
#'   sessions) the rendered report is opened in a viewer/browser. Forwarded to
#'   [quarto::quarto_render()] via its `open` handling where supported; on
#'   systems without that argument the file path is still returned.
#' @param estimands Optional character vector of estimands (parameter names) to
#'   restrict the report to. Currently informational; the template auto-detects
#'   estimands from the summary when `NULL` (default).
#'
#' @return The path to the rendered HTML file (invisibly).
#' @export
#' @seealso [performance_measures()], [plot_recovery()], [plot_coverage()],
#'   [plot_rank_ecdf()].
#' @examples
#' \dontrun{
#' result <- run_simulation(config, progress = FALSE)
#' report(result, output_file = "my-study.html")
#' }
report <- function(
  result,
  output_file = "bayesim-report.html",
  open = interactive(),
  estimands = NULL
) {
  if (!inherits(result, "bayesim_simulation_result")) {
    stop(bayesim_config_error(
      "report() requires a bayesim_simulation_result object"
    ))
  }

  rlang::check_installed("quarto", "to render simulation reports")
  if (is.null(quarto::quarto_path())) {
    stop(bayesim_config_error(paste(
      "The Quarto CLI was not found. Install it from https://quarto.org",
      "and ensure it is on your PATH, then call report() again."
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
