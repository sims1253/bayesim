# Golden workflow 1: a complete analytic simulation study.
#
# Two data conditions, a held-out test set consumed by a prediction metric,
# LinearRegressionFitter (exact conjugate posterior, no Stan), and the full
# analysis surface: performance measures, summaries, recovery/coverage plots,
# and metric column selectors. Small n and few reps keep it inside the fast
# PR gate; the assertions guard shapes and finiteness, not specific values.
library(bayesim)

.golden_gen <- function(data_spec, task_ctx) {
  n <- data_spec$n
  x <- stats::rnorm(n)
  y <- 1.5 + 0.8 * x + stats::rnorm(n, sd = 0.7)
  # Held-out evaluation rows: the pred_* metrics compare against the test
  # response and must never fall back to the training set (E2).
  x_test <- stats::rnorm(20)
  y_test <- 1.5 + 0.8 * x_test + stats::rnorm(20, sd = 0.7)
  list(
    train = data.frame(y = y, x = x),
    test = data.frame(y = y_test, x = x_test),
    response = "y",
    true_params = c(Intercept = 1.5, x = 0.8, sigma = 0.7),
    vars_of_interest = c("Intercept", "x", "sigma"),
    meta = list()
  )
}

.golden_config <- function(seed = 101L) {
  simulation_config(
    data_grid = data.frame(n = c(40L, 80L)),
    fit_grid = data.frame(model = "lm"),
    data_generator = .golden_gen,
    fitter = LinearRegressionFitter(n_draws = 200L),
    metrics = list(
      posterior_summary_metric(),
      coverage_metric(),
      pred_rmse_metric()
    ),
    n_replicates = 20L,
    seed = seed
  )
}

describe("golden workflow: complete analytic study", {
  it("runs every task, both conditions, against a held-out test set", {
    result <- run_simulation(.golden_config(), progress = FALSE)

    expect_s3_class(result, "bayesim_simulation_result")
    # 2 data conditions x 1 fit x 20 replicates, all successful.
    expect_equal(nrow(result$summary), 40L)
    expect_equal(sum(result$summary$status == "success"), 40L)
    expect_equal(length(unique(result$summary$data_n)), 2L)
  })

  it("produces finite estimator-performance measures for the main estimand", {
    result <- run_simulation(.golden_config(), progress = FALSE)
    pm <- performance_measures(result, estimand = "x")

    expect_s3_class(pm, "data.frame")
    # Fixed-truth study: the Morris, White & Crowther (2019) measures.
    measures <- c("bias", "emp_se", "mse", "coverage", "model_se", "n_sim")
    expect_true(all(measures %in% pm$measure))
    for (msr in c("bias", "emp_se", "mse", "coverage")) {
      sel <- pm$measure == msr
      expect_equal(sum(sel), 2L) # one row per data condition
      expect_true(all(is.finite(pm$value[sel])), info = msr)
      expect_false(anyNA(pm$value[sel]), info = msr)
      expect_true(all(pm$n_sim[sel] == 20L), info = msr)
    }
    expect_true(all(pm$truth_mode == "fixed"))
    # A 95% interval for a conjugate posterior should be roughly calibrated.
    sel <- pm$measure == "coverage"
    expect_true(all(pm$value[sel] > 0.7 & pm$value[sel] <= 1))
  })

  it("computes finite test-set prediction metrics", {
    result <- run_simulation(.golden_config(), progress = FALSE)
    summ <- result$summary

    expect_true(all(is.finite(summ$rmse__value)))
    # The held-out set has 20 rows; the metric must see exactly those.
    expect_true(all(summ$rmse__n_obs == 20L))
    # Plausible scale: sigma = 0.7 keeps RMSE well under 2 for n >= 40.
    expect_true(all(summ$rmse__value < 2))
  })

  it("summarizes per condition with MCSEs and sane replicate counts", {
    result <- run_simulation(.golden_config(), progress = FALSE)
    agg <- summarize_simulation(
      result,
      by = "data_n",
      metrics = metric_cols(result, "posterior_summary", fields = "mean")
    )

    expect_equal(nrow(agg), 2L)
    expect_equal(agg$n_reps, c(20L, 20L))
    expect_equal(agg$n_failed, c(0L, 0L))
    expect_equal(agg$failure_rate, c(0, 0))
    # Posterior mean of x concentrates near the truth 0.8 with shrinking MCSE
    # as n grows.
    x_mean <- paste0("posterior_summary__mean__x_", c("mean", "mcse"))
    expect_true(all(is.finite(agg[[x_mean[1]]])))
    expect_true(all(is.finite(agg[[x_mean[2]]])))
    expect_true(all(abs(agg[[x_mean[1]]] - 0.8) < 0.2))
    expect_true(agg[[x_mean[2]]][2] < agg[[x_mean[2]]][1])
  })

  it("supports recovery and coverage plots (ggplot2)", {
    skip_if_not_installed("ggplot2")
    result <- run_simulation(.golden_config(), progress = FALSE)

    expect_s3_class(plot_recovery(result, estimand = "x"), "ggplot")
    expect_s3_class(plot_coverage(result), "ggplot")
  })

  it("exposes flattened metric columns through metric_cols()", {
    result <- run_simulation(.golden_config(), progress = FALSE)

    all_cols <- metric_cols(result, "posterior_summary")
    expect_true(length(all_cols) >= 15L) # 5 fields x 3 parameters
    expect_true(all(grepl("^posterior_summary__", all_cols)))
    expect_true("posterior_summary__mean__x" %in% all_cols)

    # Field selection narrows the selector as documented.
    mean_cols <- metric_cols(result, "posterior_summary", fields = "mean")
    expect_equal(
      sort(unname(mean_cols)),
      sort(paste(
        "posterior_summary",
        "mean",
        c("Intercept", "x", "sigma"),
        sep = "__"
      ))
    )

    cov_cols <- metric_cols(result, "coverage")
    expect_true(all(grepl("^coverage__", cov_cols)))
    expect_true("coverage__by_param__x" %in% cov_cols)
  })
})
