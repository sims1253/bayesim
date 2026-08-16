# Tests for the M6 analysis & reporting layer.

# Build a synthetic summary tibble mimicking run_simulation() output.
make_summary <- function(n_per = 20, n_cond = 2) {
  rows <- list()
  for (c in seq_len(n_cond)) {
    for (r in seq_len(n_per)) {
      rows[[length(rows) + 1L]] <- data.frame(
        task_id = sprintf("task_%03d", (c - 1) * n_per + r),
        model = paste0("m", c),
        rep_idx = r,
        status = if (r == n_per) "failed" else "success",
        rmse__value = abs(rnorm(1, mean = c * 0.5, sd = 0.1)),
        coverage__mean = c(1, 0)[1 + (r %% 2)],
        rank__n_draws = 100L,
        rank__by_param__b_x = as.integer(runif(1, 0, 100)),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

describe("summarize_simulation", {
  it("uses field roles to exclude counts from scientific aggregates", {
    df <- data.frame(
      model = "x",
      status = c("success", "success"),
      pred__value = c(1, 2),
      pred__n_obs = c(10, 20)
    )
    agg <- summarize_simulation(
      df,
      by = "model",
      metrics = c("pred__value", "pred__n_obs")
    )
    expect_true("pred__value_mean" %in% names(agg))
    expect_false("pred__n_obs_mean" %in% names(agg))
  })

  it("rejects unknown grouping and metric columns", {
    df <- data.frame(model = "x", value__x = 1)
    expect_error(
      summarize_simulation(df, by = "missing"),
      class = "bayesim_config_error"
    )
    expect_error(
      summarize_simulation(df, metrics = "missing__value"),
      class = "bayesim_config_error"
    )
  })

  it("aggregates metrics with mean/median/sd/mcse per condition", {
    df <- make_summary(n_per = 20, n_cond = 2)
    agg <- summarize_simulation(df, by = "model", metrics = "rmse__value")
    expect_s3_class(agg, "data.frame")
    expect_equal(nrow(agg), 2L)
    expect_true("rmse__value_mean" %in% names(agg))
    expect_true("rmse__value_mcse" %in% names(agg))
    expect_true("n_reps" %in% names(agg))
    # MCSE is finite and positive for the rmse column.
    expect_true(all(is.finite(agg$rmse__value_mcse)))
  })

  it("computes failure rate from status", {
    df <- make_summary(n_per = 20, n_cond = 1)
    agg <- summarize_simulation(df, by = "model")
    expect_equal(agg$n_failed, 1L)
    expect_equal(agg$failure_rate, 1 / 20)
  })

  it("does not treat task and metric errors as condition columns", {
    df <- data.frame(
      task_id = c("t1", "t2", "t3", "t4"),
      model = rep(c("a", "b"), each = 2L),
      status = c("success", "failed", "success", "success"),
      error_class = c(NA, "bayesim_fit_error", NA, NA),
      error_message = c(NA, "fit failed", NA, NA),
      optional__error_class = c(NA, "simpleError", NA, NA),
      optional__error_message = c(NA, "metric failed", NA, NA),
      rmse__value = c(1, NA, 2, 4),
      stringsAsFactors = FALSE
    )

    agg <- summarize_simulation(df)

    expect_equal(nrow(agg), 2L)
    expect_equal(agg$n_reps[agg$model == "a"], 2L)
    expect_equal(agg$n_failed[agg$model == "a"], 1L)
    expect_equal(agg$failure_rate[agg$model == "a"], 0.5)
    expect_false(any(grepl("error_(class|message)", names(agg))))
  })

  it("reports how many finite values contributed to each metric", {
    df <- data.frame(
      model = "x",
      status = "success",
      estimate__value = c(1, 2, NA, Inf)
    )

    agg <- summarize_simulation(df, by = "model", metrics = "estimate__value")

    expect_equal(agg$estimate__value_n_used, 2L)
    expect_equal(agg$estimate__value_mean, 1.5)
  })

  it("coverage MCSE uses sqrt(p(1-p)/n)", {
    # coverage values all 1 -> p=1 -> MCSE = sqrt(1*0/n) = 0
    df <- data.frame(
      model = "x",
      status = "success",
      coverage__mean = rep(1, 10)
    )
    agg <- summarize_simulation(df, by = "model", metrics = "coverage__mean")
    expect_equal(agg$coverage__mean_mcse, 0)
  })

  it("honors field-level schema metadata when aggregating a result", {
    df <- data.frame(
      model = "x",
      status = c("success", "success"),
      custom__estimate = c(1, 3),
      custom__count = c(10, 20),
      check.names = FALSE
    )
    result <- structure(
      list(
        summary = df,
        metric_field_metadata = list(
          custom = list(
            estimate = list(
              role = "estimate",
              aggregation = "mean",
              mcse = "sd"
            ),
            count = list(role = "count", aggregation = "none", mcse = "none")
          )
        )
      ),
      class = "bayesim_simulation_result"
    )
    agg <- summarize_simulation(result, by = "model")
    expect_true("custom__estimate_mean" %in% names(agg))
    expect_false("custom__count_mean" %in% names(agg))
  })

  it("coverage MCSE for p=0.5, n=100 is sqrt(0.25/100) = 0.05", {
    df <- data.frame(
      model = "x",
      status = "success",
      coverage__mean = rep(c(0, 1), 50)
    )
    agg <- summarize_simulation(df, by = "model", metrics = "coverage__mean")
    expect_equal(agg$coverage__mean_mcse, sqrt(0.25 / 100))
  })

  it("bias MCSE = sd/sqrt(n) for general metrics", {
    vals <- c(1, 2, 3, 4)
    df <- data.frame(model = "x", status = "success", rmse__value = vals)
    agg <- summarize_simulation(df, by = "model", metrics = "rmse__value")
    # rmse-like -> delta-method, but verify the general bias path separately.
    # Here test a non-rmse metric column name.
    df2 <- data.frame(model = "x", status = "success", estimate__value = vals)
    agg2 <- summarize_simulation(df2, by = "model", metrics = "estimate__value")
    expect_equal(agg2$estimate__value_mcse, stats::sd(vals) / sqrt(4))
  })

  it("uses collision-safe group identities", {
    df <- data.frame(
      left = c("a|b", "a"),
      right = c("c", "b|c"),
      value__x = c(1, 2)
    )
    expect_equal(nrow(summarize_simulation(df, c("left", "right"))), 2L)

    na_df <- data.frame(group = c(NA, "NA"), value__x = c(1, 2))
    expect_equal(nrow(summarize_simulation(na_df, "group")), 2L)
  })

  it("groups factors and mixed numeric/character columns", {
    df <- data.frame(
      fac = factor(c("a", "a", "b")),
      number = c(1, 1, 2),
      label = c("1", "1", "2"),
      value__x = c(1, 3, 5)
    )
    agg <- summarize_simulation(df, c("fac", "number", "label"))
    expect_equal(nrow(agg), 2L)
    expect_equal(agg$value__x_mean, c(2, 5))
  })

  it("stays silent programmatically even for 100+ metric columns", {
    n_cols <- 120L
    wide <- data.frame(
      model = rep(c("a", "b"), each = 5L),
      status = "success"
    )
    for (i in seq_len(n_cols)) {
      wide[[sprintf("m%03d__value", i)]] <- i + seq_len(10) * 0.1
    }
    agg <- expect_silent(summarize_simulation(wide, by = "model"))
    expect_equal(nrow(agg), 2L)
    # Nothing is dropped or truncated to make the summary discoverable.
    expect_equal(sum(grepl("_mean$", names(agg))), n_cols)
  })
})

describe("wide summary discoverability hint", {
  it("fires only interactively, only when wide, only by default", {
    # Interactive + wide + no explicit metrics -> hint mentions narrowing.
    expect_message(
      maybe_wide_summary_hint(30L, NULL, .interactive = TRUE),
      "metrics"
    )
    # Noninteractive programmatic use is completely silent, however wide.
    expect_silent(maybe_wide_summary_hint(30L, NULL, .interactive = FALSE))
    # User already narrowed via metrics = -> no hint even interactively.
    expect_silent(
      maybe_wide_summary_hint(30L, "rmse__value", .interactive = TRUE)
    )
    # Narrow summaries never hint.
    expect_silent(maybe_wide_summary_hint(3L, NULL, .interactive = TRUE))
  })

  it("threshold matches the documented wide-summary cutoff", {
    expect_silent(
      maybe_wide_summary_hint(
        WIDE_SUMMARY_HINT_THRESHOLD,
        NULL,
        .interactive = TRUE
      )
    )
    expect_message(
      maybe_wide_summary_hint(
        WIDE_SUMMARY_HINT_THRESHOLD + 1L,
        NULL,
        .interactive = TRUE
      )
    )
  })
})

describe("metric_cols", {
  df <- data.frame(
    task_id = c("t1", "t2"),
    posterior_summary__mean__x = c(1, 2),
    posterior_summary__sd__x = c(0.1, 0.2),
    posterior_summary__n_eff = c(100, 110),
    coverage__value__x = c(TRUE, FALSE),
    check.names = FALSE
  )

  it("returns named flattened column names and filters fields", {
    cols <- metric_cols(df, "posterior_summary", fields = "mean")
    expect_equal(
      cols,
      c(mean__x = "posterior_summary__mean__x")
    )
  })

  it("returns long metric values with optional parameter", {
    long <- metric_cols(df, "posterior_summary", as = "long")
    expect_named(long, c("task_id", "field", "param", "value"))
    expect_equal(nrow(long), 6L)
    expect_true(all(is.na(long$param[long$field == "n_eff"])))
    expect_equal(long$param[long$field == "mean"], rep("x", 2L))
  })

  it("accepts a result object and reports available prefixes", {
    result <- structure(list(summary = df), class = "bayesim_simulation_result")
    expect_length(metric_cols(result, "coverage"), 1L)
    expect_error(
      metric_cols(df, "missing"),
      "coverage, posterior_summary",
      class = "bayesim_config_error"
    )
  })
})

describe("plot_recovery estimand/var terminology", {
  recovery_result <- function() {
    summary <- data.frame(
      status = rep("success", 4),
      truth__x = rep(1, 4),
      posterior_summary__mean__x = c(0.9, 1.1, 1.0, 1.2),
      posterior_summary__sd__x = rep(0.2, 4),
      posterior_summary__q_lower__x = rep(0.6, 4),
      posterior_summary__q_upper__x = rep(1.4, 4),
      check.names = FALSE
    )
    structure(
      list(summary = summary),
      class = "bayesim_simulation_result"
    )
  }

  it("accepts estimand as the primary argument name", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    p <- plot_recovery(recovery_result(), estimand = "x")
    expect_s3_class(p, "ggplot")
    expect_equal(p$labels$x, "true x")
    expect_equal(p$labels$y, "posterior mean x")
  })

  it("keeps the legacy var alias and produces the same plot", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    p_var <- plot_recovery(recovery_result(), var = "x")
    p_estimand <- plot_recovery(recovery_result(), estimand = "x")
    expect_s3_class(p_var, "ggplot")
    expect_equal(p_var$labels, p_estimand$labels)
    expect_equal(
      vapply(p_var$layers, function(l) class(l$geom)[1], character(1)),
      vapply(p_estimand$layers, function(l) class(l$geom)[1], character(1))
    )
  })

  it("errors clearly when neither estimand nor var is provided", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    expect_error(
      plot_recovery(recovery_result()),
      "estimand",
      class = "bayesim_config_error"
    )
  })

  it("errors when estimand and var disagree", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    expect_error(
      plot_recovery(recovery_result(), estimand = "x", var = "y"),
      class = "bayesim_config_error"
    )
  })

  it("accepts identical estimand and var without error", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    p <- plot_recovery(recovery_result(), estimand = "x", var = "x")
    expect_s3_class(p, "ggplot")
  })
})

describe("sbc_ranks", {
  it("extracts ranks into a long tibble", {
    df <- make_summary(n_per = 5, n_cond = 1)
    ranks <- sbc_ranks(df)
    expect_s3_class(ranks, "data.frame")
    expect_true(all(c("task_id", "param", "rank", "n_draws") %in% names(ranks)))
    expect_true(all(ranks$param == "b_x"))
    expect_equal(nrow(ranks), 5L)
  })

  it("surfaces n_ranks per variable when present (F4)", {
    df <- make_summary(n_per = 4, n_cond = 1)
    # Add the F4 per-variable n_ranks column.
    df[["rank__n_ranks__b_x"]] <- 101L
    ranks <- sbc_ranks(df)
    expect_true("n_ranks" %in% names(ranks))
    expect_true(all(ranks$n_ranks == 101L))
  })

  it("preserves simulation condition columns for faceting", {
    df <- make_summary(n_per = 4, n_cond = 2)
    df$data_n <- rep(c(20L, 40L), each = 4L)

    ranks <- sbc_ranks(df)

    expect_equal(unique(ranks$data_n), c(20L, 40L))
  })

  it("returns empty tibble when no rank columns", {
    df <- data.frame(task_id = "t1", model = "x", status = "success")
    ranks <- sbc_ranks(df)
    expect_equal(nrow(ranks), 0L)
  })

  it("handles single-param flattened rank column (no __param suffix)", {
    # The flatten layer collapses length-1 named vectors to a scalar
    # `rank__by_param`; sbc_ranks must still extract it.
    df <- data.frame(
      task_id = c("t1", "t2"),
      model = "x",
      status = "success",
      rank__n_draws = c(100L, 100L),
      rank__by_param = c(40L, 60L),
      check.names = FALSE
    )
    ranks <- sbc_ranks(df)
    expect_equal(nrow(ranks), 2L)
    expect_true(all(c("task_id", "param", "rank") %in% names(ranks)))
  })
})

describe("plot functions", {
  it("plot_coverage uses the nominal level of the intervals it summarizes", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    summary <- data.frame(
      status = rep("success", 4),
      truth__x = rep(0, 4),
      posterior_summary__mean__x = c(-0.1, 0, 0.1, 0),
      posterior_summary__q_lower__x = rep(-1, 4),
      posterior_summary__q_upper__x = rep(1, 4),
      coverage__by_param__x = rep(1, 4),
      check.names = FALSE
    )
    result <- structure(
      list(
        summary = summary,
        metric_field_metadata = list(
          posterior_summary = list(
            q_lower = list(nominal = 0.95),
            q_upper = list(nominal = 0.95)
          ),
          coverage = list(by_param = list(nominal = 0.90))
        )
      ),
      class = "bayesim_simulation_result"
    )
    built <- ggplot2::ggplot_build(plot_coverage(result))
    expect_equal(unique(built$data[[1]]$yintercept), 0.95)
  })

  it("plot_rank_hist returns a ggplot object", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    df <- make_summary(n_per = 10, n_cond = 1)
    p <- plot_rank_hist(sbc_ranks(df))
    expect_s3_class(p, "ggplot")
  })

  it("plot_rank_ecdf returns a ggplot with ribbon + line", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    df <- make_summary(n_per = 10, n_cond = 1)
    p <- plot_rank_ecdf(sbc_ranks(df))
    expect_s3_class(p, "ggplot")
    # The ribbon geom should be present.
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true("GeomRibbon" %in% geoms)
  })

  it("plot_rank_ecdf computes separate condition facets", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    df <- make_summary(n_per = 10, n_cond = 2)
    df$data_n <- rep(c(20L, 40L), each = 10L)

    p <- plot_rank_ecdf(sbc_ranks(df), by = "data_n")

    expect_s3_class(p, "ggplot")
    expect_equal(sort(unique(p$data$data_n)), c(20L, 40L))
  })

  it("plot_metric returns a ggplot", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    df <- make_summary(n_per = 10, n_cond = 1)
    p <- plot_metric(df, "rmse__value", x = "task_id")
    expect_s3_class(p, "ggplot")
  })

  it("plot_metric errors on unknown metric", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    df <- make_summary(n_per = 5, n_cond = 1)
    expect_error(plot_metric(df, "nonexistent"), class = "bayesim_config_error")
  })

  it("plot_metric validates requested x and facet columns", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    df <- make_summary(n_per = 5, n_cond = 1)
    expect_error(
      plot_metric(df, "rmse__value", x = "missing"),
      class = "bayesim_config_error"
    )
    expect_error(
      plot_metric(df, "rmse__value", facets = "missing"),
      class = "bayesim_config_error"
    )
  })
})

describe("plot_rank_ecdf per-task rank normalization", {
  it("normalizes each task by its own support when supports are mixed", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    set.seed(20260816)
    n_small <- 50L
    n_big <- 50L
    ranks <- tibble::tibble(
      task_id = sprintf("t%03d", seq_len(n_small + n_big)),
      param = "x",
      rank = as.integer(c(
        sample(0:100, n_small, replace = TRUE),
        sample(0:200, n_big, replace = TRUE)
      )),
      n_draws = 100L,
      n_ranks = as.integer(c(rep(101L, n_small), rep(201L, n_big)))
    )

    # Mixed finite n_ranks within a panel -> one warning about per-task
    # normalization and the approximate band.
    expect_warning(plot_rank_ecdf(ranks), "normalized per task")
    p <- suppressWarnings(plot_rank_ecdf(ranks))

    d <- p$data
    # Normalized values live strictly inside (0, 1).
    expect_true(all(d$rank_norm > 0 & d$rank_norm < 1))
    # Uniform ranks on each task's own support -> the ECDF at 0.5 stays near
    # 0.5. Under the old group-max normalization the small-support tasks
    # squash toward 0 and the ECDF at 0.5 sits near 0.75.
    expect_lt(abs(mean(d$rank_norm <= 0.5) - 0.5), 0.15)
    # Exact check: every row was scaled by its own support.
    expect_equal(d$rank_norm, sort((ranks$rank + 0.5) / ranks$n_ranks))
    # Small-support tasks reach near 1, far beyond the (S_small + 0.5) /
    # S_big ceiling that a shared group-max normalization would impose.
    small_norm <- (ranks$rank[ranks$n_ranks == 101L] + 0.5) / 101
    expect_gt(max(small_norm), (100 + 0.5) / 200)
  })

  it("uses exact per-row normalized values in a deterministic mixed case", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    ranks <- tibble::tibble(
      task_id = c("a", "b", "c"),
      param = "x",
      rank = c(0L, 50L, 199L),
      n_draws = c(100L, 100L, 200L),
      n_ranks = c(101L, 101L, 201L)
    )

    expect_warning(plot_rank_ecdf(ranks), "normalized per task")
    p <- suppressWarnings(plot_rank_ecdf(ranks))

    d <- p$data
    expect_equal(
      d$rank_norm,
      sort(c(0.5 / 101, 50.5 / 101, 199.5 / 201))
    )
    expect_equal(d$ecdf, seq_len(3L) / 3L)
    # The band grid uses the largest support in the panel.
    expect_equal(d$S, rep(200, 3L))
  })

  it("matches the golden pooled formula when all tasks share one support", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    ranks <- tibble::tibble(
      task_id = sprintf("t%d", 1:6),
      param = "x",
      rank = c(3L, 17L, 42L, 8L, 90L, 55L),
      n_draws = 100L,
      n_ranks = 101L
    )

    # Homogeneous support: no warning, byte-identical pooled normalization.
    p <- expect_silent(plot_rank_ecdf(ranks))

    d <- p$data
    golden <- sort((c(3, 17, 42, 8, 90, 55) + 0.5) / 101)
    expect_equal(d$rank_norm, golden)
    expect_equal(d$ecdf, seq_len(6L) / 6L)
    expect_equal(d$S, rep(100, 6L))
  })

  it("falls back per row to n_draws for legacy results without n_ranks", {
    skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
    ranks <- tibble::tibble(
      task_id = c("a", "b", "c"),
      param = "x",
      rank = c(10L, 150L, 50L),
      n_draws = c(100L, 200L, NA_integer_),
      n_ranks = c(NA_integer_, NA_integer_, NA_integer_)
    )

    # No finite n_ranks anywhere: no heterogeneity warning; each row uses its
    # own n_draws (+1), and the row without even n_draws falls back to the
    # panel max rank.
    p <- expect_silent(plot_rank_ecdf(ranks))

    expect_equal(
      p$data$rank_norm,
      sort(c(10.5 / 101, 150.5 / 201, 50.5 / 151))
    )
  })
})

describe("SBC simultaneous bands", {
  it("rejects non-probabilities and invalid dimensions", {
    expect_error(sbc_band(10, conf_level = 0), class = "bayesim_config_error")
    expect_error(adjust_gamma(10, L = 1, K = 0), class = "bayesim_config_error")
  })
  it("returns a valid envelope with fixed endpoints", {
    band <- sbc_band(N = 25L, K = 10L, conf_level = 0.95)

    expect_length(band$x, 11L)
    expect_equal(band$x[c(1L, 11L)], c(0, 1))
    expect_true(all(band$lower <= band$upper))
    expect_true(all(diff(band$lower) >= 0))
    expect_true(all(diff(band$upper) >= 0))
  })

  it("handles the one-interval degenerate band", {
    band <- sbc_band(N = 10L, K = 1L, conf_level = 0.95)
    expect_length(band$x, 2L)
    expect_equal(band$lower[c(1L, 2L)], c(0, 1))
    expect_equal(band$upper[c(1L, 2L)], c(0, 1))
  })

  it("matches the ported reference gamma for a fixed design", {
    expect_equal(
      adjust_gamma(20L, L = 1L, K = 20L, conf_level = 0.95),
      0.0140491762611022,
      tolerance = 1e-12
    )
  })

  it("announces the conservative fallback for multiple chains", {
    # warn_once(): the fallback is announced as a warning (once per run)
    # rather than a message on every call.
    expect_warning(
      multi <- adjust_gamma(20L, L = 2L, K = 20L, conf_level = 0.95),
      "single-sample SBC band"
    )
    single <- adjust_gamma(20L, L = 1L, K = 20L, conf_level = 0.95)
    expect_equal(multi, single)
  })
})
