# E3: performance_measures() — Morris et al. (2019) estimator-performance layer.
skip_on_cran()

.gen <- function(data_spec, task_ctx) {
  n <- data_spec$n; b <- data_spec$beta
  x <- stats::rnorm(n)
  y <- 1 + b * x + stats::rnorm(n)
  list(train = data.frame(y = y, x = x), test = NULL, response = "y",
       true_params = c(Intercept = 1, x = b, sigma = 1),
       vars_of_interest = c("Intercept", "x", "sigma"),
       meta = list())
}

describe("performance_measures", {
  it("returns a tidy tibble of measures with MCSE per estimand/condition", {
    config <- simulation_config(
      data_grid = data.frame(n = 60L, beta = c(0.5, 1.0)),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 300L),
      metrics = list(posterior_summary_metric(), coverage_metric()),
      n_replicates = 30L, seed = 7L
    )
    res <- run_simulation(config, resume = "never", progress = FALSE)
    pm <- performance_measures(res)

    expect_s3_class(pm, "data.frame")
    expect_true(all(c("estimand", "measure", "value", "mcse", "n_sim") %in% names(pm)))
    measures <- unique(pm$measure)
    expect_true("bias" %in% measures)
    expect_true("emp_se" %in% measures)
    expect_true("mse" %in% measures)
    expect_true("coverage" %in% measures)
    expect_true("model_se" %in% measures)
    expect_true("n_sim" %in% measures)
    # All n_sim rows equal the replicate count (30) per condition.
    nsim <- pm[pm$measure == "n_sim", ]
    expect_true(all(nsim$n_sim == 30L))
  })

  it("coverage is near nominal (0.9) for the 'x' coefficient", {
    config <- simulation_config(
      data_grid = data.frame(n = 100L, beta = 0.5),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 500L),
      metrics = list(posterior_summary_metric()),
      n_replicates = 100L, seed = 99L
    )
    res <- run_simulation(config, resume = "never", progress = FALSE)
    pm <- performance_measures(res, estimand = "x")
    cov <- pm$value[pm$measure == "coverage"]
    # 95% interval; generous band for MC noise (nominal 0.95).
    expect_true(cov > 0.85 && cov < 1.0)
  })

  it("supports a single estimand and a custom point estimator", {
    config <- simulation_config(
      data_grid = data.frame(n = 60L, beta = 0.5),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 200L),
      metrics = list(posterior_summary_metric()),
      n_replicates = 10L, seed = 3L
    )
    res <- run_simulation(config, resume = "never", progress = FALSE)
    pm <- performance_measures(res, estimand = "x", estimator = "median")
    expect_true(all(pm$estimand == "x"))
    expect_true("bias" %in% pm$measure)
  })

  it("errors when no truth/posterior-summary columns are present", {
    df <- data.frame(task_id = "a", rep_idx = 1L, status = "success",
                     posterior_summary__mean__x = 0.5, fit_model = "m")
    expect_error(performance_measures(df), class = "bayesim_config_error")
  })
})
