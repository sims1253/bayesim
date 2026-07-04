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
#'   MCSE formulas follow rsimsum (Gasparini, 2018):
#'   \itemize{
#'     \item bias / mean metrics: MCSE = sd / sqrt(n)
#'     \item coverage (0/1) metrics: MCSE = sqrt(p(1-p) / n)
#'     \item rmse metrics: MCSE via the delta method on the mean and variance of
#'       the squared errors (see Details).
#'   }
#'
#' @param result A `bayesim_simulation_result` object (uses `result$summary`),
#'   or a data.frame/tibble of per-task metrics.
#' @param by Character vector of grouping columns (conditions). Defaults to the
#'   data_grid and fit_grid columns found in the summary (excluding `task_id`,
#'   `rep_idx`, `status`, and metric columns).
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

  meta_cols <- c("task_id", "rep_idx", "status", "timing_total", "timing_warmup", "timing_sample")
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  # exclude timing/pure-metadata numeric columns from default metrics
  exclude_meta <- numeric_cols[grepl("^(timing|n_)|(_timing)$", numeric_cols)]
  metric_cols <- setdiff(numeric_cols, c(meta_cols, exclude_meta, by %||% character()))

  if (!is.null(metrics)) metric_cols <- intersect(metric_cols, metrics)
  if (length(metric_cols) == 0L) {
    warning("No metric columns found to summarize")
    return(tibble::as_tibble(df[0, , drop = FALSE]))
  }

  if (is.null(by)) {
    # Auto-detect grouping columns: everything non-numeric that isn't metadata,
    # plus any explicit condition columns.
    candidate <- names(df)[!vapply(df, is.numeric, logical(1))]
    by <- setdiff(candidate, c("task_id", "status"))
  }
  by <- intersect(by, names(df))

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
      out[[paste0(m, "_median")]] <- if (n > 0L) stats::median(vals) else NA_real_
      sd_v <- if (n > 1L) stats::sd(vals) else NA_real_
      out[[paste0(m, "_sd")]] <- sd_v
      # E4: MCSE by summary_type. Coverage-like columns (binary 0/1 or a column
      # the metric declares as "proportion") use sqrt(p(1-p)/n). Everything
      # else (including pred-RMSE means) uses the plain sd/sqrt(n) MCSE for what
      # is reported (the mean of per-task values) — the previous delta-method
      # rmse branch re-squared per-task RMSEs, producing an MCSE for a quantity
      # (sqrt-mean-MSE) that is not the reported _mean column.
      if (n > 1L) {
        if (is_coverage_like(m, vals)) {
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

#' Heuristic: is a metric column a coverage rate (binary or between 0 and 1)?
#' @keywords internal
is_coverage_like <- function(name, vals) {
  grepl("coverage", name, ignore.case = TRUE) ||
    (all(vals >= 0 & vals <= 1, na.rm = TRUE) && length(unique(vals)) <= 3L)
}

#' Heuristic: is a metric column an RMSE-like quantity (squared-error-based)?
#' @keywords internal
is_rmse_like <- function(name) {
  grepl("rmse|mse", name, ignore.case = TRUE)
}

#' Delta-method MCSE for RMSE = sqrt(mean(e^2)).
#' Var(sqrt(M)) approximates Var(e^2) / (4*M) where M = mean(e^2).
#' MCSE = sqrt(Var(e^2) / (4 * n * M)).
#' @keywords internal
rmse_mcse <- function(vals) {
  n <- length(vals)
  if (n < 2L) return(NA_real_)
  e2 <- vals^2
  M <- mean(e2)
  if (M <= 0) return(NA_real_)
  sqrt(stats::var(e2) / (4 * n * M))
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
    task_id = character(0), param = character(0),
    rank = integer(0), n_draws = integer(0), n_ranks = integer(0)
  )
  if (is.null(df)) return(empty)
  # rank__by_param__<param> columns hold the per-task per-parameter rank counts.
  # For the single-parameter case the flattener collapses the length-1 named
  # vector to a scalar `rank__by_param` (no __<param> suffix); handle both.
  rank_cols_multi <- grep("^rank__by_param__", names(df), value = TRUE)
  rank_cols_single <- if ("rank__by_param" %in% names(df)) "rank__by_param" else character(0)
  rank_cols <- c(rank_cols_multi, rank_cols_single)
  n_draws_col <- grep("^rank__n_draws$", names(df), value = TRUE)
  if (length(rank_cols) == 0L) return(empty)
  n_draws <- if (length(n_draws_col) == 1L) df[[n_draws_col]] else NA_integer_
  rows <- list()
  for (col in rank_cols) {
    param <- sub("^rank__by_param__?", "", col)
    if (param == "") param <- "(single)"
    # Per-variable n_ranks when present (F4); fall back to n_draws + 1.
    n_ranks_col <- paste0("rank__n_ranks__", param)
    n_ranks <- if (n_ranks_col %in% names(df)) as.integer(df[[n_ranks_col]]) else NA_integer_
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
  if (inherits(ranks, "bayesim_simulation_result")) ranks <- sbc_ranks(ranks)
  if (nrow(ranks) == 0L) stop(bayesim_config_error("No rank data; was rank_metric() used?"))
  ggplot2::ggplot(ranks, ggplot2::aes(.data$rank)) +
    ggplot2::geom_histogram(bins = 20, color = "white") +
    ggplot2::facet_wrap(~param, scales = "free_x") +
    ggplot2::labs(x = "rank", y = "count", title = "SBC rank histograms") +
    ggplot2::theme_minimal()
}

#' Plot SBC rank ECDF with simultaneous confidence bands
#'
#' @description Plots the empirical CDF of SBC ranks against the uniform CDF,
#'   with simultaneous confidence bands following Säilynoja, Bürkner & Vehtari
#'   (2022). Under correct calibration the ECDF should stay within the band.
#' @param ranks A tibble from [sbc_ranks()], or a `bayesim_simulation_result`.
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' plot_rank_ecdf(sbc_ranks(result))
#' }
plot_rank_ecdf <- function(ranks) {
  rlang::check_installed("ggplot2", "to use plot_rank_ecdf()")
  if (inherits(ranks, "bayesim_simulation_result")) ranks <- sbc_ranks(ranks)
  if (nrow(ranks) == 0L) stop(bayesim_config_error("No rank data; was rank_metric() used?"))

  # Build ECDF and simultaneous confidence band per parameter (Säilynoja,
  # Bürkner & Vehtari 2022), ported via adjust_gamma()/sbc_band().
  conf_level <- 0.95
  plot_data <- do.call(rbind, lapply(unique(ranks$param), function(p) {
    sub <- ranks[ranks$param == p, , drop = FALSE]
    n <- nrow(sub)
    # Prefer post-thinning n_ranks (F4); fall back to n_draws for old results.
    S <- NA
    if ("n_ranks" %in% names(sub) && any(is.finite(sub$n_ranks))) {
      S <- max(sub$n_ranks - 1L, na.rm = TRUE)
    }
    if (!is.finite(S) || S < 1L) S <- max(sub$n_draws, na.rm = TRUE)
    if (!is.finite(S) || S < 1L) S <- max(sub$rank, na.rm = TRUE)
    # Normalized rank on [0,1]: rank in 0..S, map to (rank+0.5)/(S+1).
    r <- (sort(sub$rank) + 0.5) / (S + 1)
    ecdf_y <- seq_len(n) / n
    # Simultaneous band over the ECDF of n ranks, using the exact single-sample
    # construction (sbc_band -> adjust_gamma_optimize).
    band <- tryCatch(
      sbc_band(N = n, K = max(100L, n), conf_level = conf_level),
      error = function(e) NULL
    )
    if (is.null(band)) {
      # Degenerate fallback (very small n): a wide pointwise band.
      bound <- sqrt(log(2 / (1 - conf_level)) / (2 * n))
      tibble::tibble(
        param = p,
        rank_norm = r,
        ecdf = ecdf_y,
        lower = pmax(0, ecdf_y - bound),
        upper = pmin(1, ecdf_y + bound)
      )
    } else {
      # Interpolate the band onto the observed rank points.
      lo <- stats::approx(band$x, band$lower, xout = r, rule = 2)$y
      hi <- stats::approx(band$x, band$upper, xout = r, rule = 2)$y
      tibble::tibble(
        param = p,
        rank_norm = r,
        ecdf = ecdf_y,
        lower = pmax(0, lo),
        upper = pmin(1, hi)
      )
    }
  }))

  ggplot2::ggplot(plot_data, ggplot2::aes(.data$rank_norm)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      alpha = 0.2, fill = "grey50"
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$ecdf)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    ggplot2::facet_wrap(~param) +
    ggplot2::labs(
      x = "normalized rank", y = "ECDF",
      title = "SBC rank ECDF with 95% simultaneous confidence band"
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
      "Posterior summary column " %+% mean_col %+% " not found; "
        %+% "compute posterior_summary_metric() first."
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
    truth = if (!is.na(truth_col) && truth_col %in% names(df)) df[[truth_col]] else NA_real_,
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
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::labs(
      x = paste0("true ", var), y = paste0("posterior mean ", var),
      title = paste0("Parameter recovery: ", var)
    ) +
    ggplot2::theme_minimal()
  # Posterior-interval segments by default (E7).
  if (lower_col %in% names(df)) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), alpha = 0.2)
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
    ggplot2::geom_hline(yintercept = nominal, color = "grey50", linetype = "dashed") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$value - .data$mcse,
                                        ymax = .data$value + .data$mcse), width = 0.2) +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "estimand", y = paste0("coverage (nominal ", nominal, ")"),
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
    stop(bayesim_config_error("Metric column " %+% metric %+% " not found in summary"))
  }
  df$.metric <- df[[metric]]
  if (is.null(x)) x <- "task_id"
  p <- ggplot2::ggplot(df, ggplot2::aes(.data[[x]], .data$.metric)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::labs(y = metric, title = metric) +
    ggplot2::theme_minimal()
  if (!is.null(facets)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facets)))
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
performance_measures <- function(result, estimand = NULL, estimator = c("mean", "median"),
                                 by = NULL) {
  estimator <- match.arg(estimator)
  df <- if (is.data.frame(result)) result else result$summary
  if (is.null(df) || !is.data.frame(df)) {
    stop(bayesim_config_error(
      "performance_measures requires a simulation result or summary data.frame"
    ))
  }

  # Default grouping: data_grid/fit_grid condition columns. Exclude task
  # metadata, per-task diagnostics, metric columns (carry __), truth columns,
  # and timing/n_-prefixed columns.
  if (is.null(by)) {
    exclude_patterns <- "^(task_id|rep_idx|status|timing|fit_model|rhat|ess|divergent|max_treedepth|n_)|__|^truth__"
    by <- names(df)[!grepl(exclude_patterns, names(df))]
    # Keep only scalar condition columns.
    by <- by[vapply(df[, by, drop = FALSE], function(col) {
      is.numeric(col) || is.character(col) || is.factor(col) || is.logical(col)
    }, logical(1))]
  }
  if (length(by) == 0L) by <- "fit_model"

  # Discover estimands: parameters with both truth__ and posterior_summary__mean__.
  truth_cols <- grep("^truth__", names(df), value = TRUE)
  est_col_tmpl <- paste0("posterior_summary__", estimator, "__")
  available_estimands <- vapply(truth_cols, function(tc) {
    param <- sub("^truth__", "", tc)
    paste0(est_col_tmpl, param) %in% names(df)
  }, logical(1))
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

    if (!truth_col %in% names(df) || !est_col %in% names(df)) next

    split_cols <- if (length(by)) df[, by, drop = FALSE] else data.frame(.group = rep(1, nrow(df)))
    if (length(by) == 0L) by_use <- ".group" else by_use <- by
    splits <- interaction(split_cols, drop = TRUE)

    for (lev in levels(splits)) {
      sel <- which(splits == lev)
      est <- df[[est_col]][sel]
      truth <- df[[truth_col]][sel]
      ok <- !is.na(est) & !is.na(truth)
      est <- est[ok]; truth <- truth[ok]
      n <- length(est)
      if (n == 0L) next
      errs <- est - truth

      # Condition cell values for the `by` columns.
      cond <- as.list(df[sel[1L], by_use, drop = FALSE])

      add <- function(measure, value, mcse) {
        row <- c(cond, list(estimand = param, measure = measure,
                            value = value, mcse = mcse, n_sim = n))
        rows[[length(rows) + 1L]] <<- row
      }

      # bias = mean(est - truth); MCSE = sd(est)/sqrt(n)
      add("bias", mean(errs), if (n > 1L) stats::sd(est) / sqrt(n) else NA_real_)
      # empirical SE = sd(est); MCSE = sd/sqrt(2(n-1))
      emp_se <- if (n > 1L) stats::sd(est) else NA_real_
      add("emp_se", emp_se, if (n > 2L) emp_se / sqrt(2 * (n - 1L)) else NA_real_)
      # MSE = mean((est-truth)^2); MCSE = sqrt(Var(err^2)/n)
      sq <- errs^2
      add("mse", mean(sq), if (n > 1L) sqrt(stats::var(sq) / n) else NA_real_)
      # average model SE = mean(posterior_sd); MCSE = sd(posterior_sd)/sqrt(n)
      if (sd_col %in% names(df)) {
        pse <- df[[sd_col]][sel][ok]
        add("model_se", mean(pse, na.rm = TRUE),
            if (sum(!is.na(pse)) > 1L) stats::sd(pse, na.rm = TRUE) / sqrt(sum(!is.na(pse))) else NA_real_)
      }
      # coverage = P(truth in [lo, hi]); MCSE = sqrt(p(1-p)/n)
      if (lo_col %in% names(df) && hi_col %in% names(df)) {
        lo <- df[[lo_col]][sel][ok]; hi <- df[[hi_col]][sel][ok]
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

  tibble::as_tibble(do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE)))
}
