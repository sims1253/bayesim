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

  it("coverage MCSE uses sqrt(p(1-p)/n)", {
    # coverage values all 1 -> p=1 -> MCSE = sqrt(1*0/n) = 0
    df <- data.frame(
      model = "x", status = "success",
      coverage__mean = rep(1, 10)
    )
    agg <- summarize_simulation(df, by = "model", metrics = "coverage__mean")
    expect_equal(agg$coverage__mean_mcse, 0)
  })

  it("coverage MCSE for p=0.5, n=100 is sqrt(0.25/100) = 0.05", {
    df <- data.frame(
      model = "x", status = "success",
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
      model = "x", status = "success",
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
})
