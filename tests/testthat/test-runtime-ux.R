# F1-F4: runtime UX helpers.
skip_on_cran()

.gen <- function(data_spec, task_ctx) {
  n <- data_spec$n
  x <- stats::rnorm(n)
  y <- x + stats::rnorm(n)
  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(x = 1),
    vars_of_interest = "x",
    meta = list()
  )
}

describe("F1 preflight", {
  it("reports task count, grid shape, and unmet needs", {
    config <- simulation_config(
      data_grid = data.frame(n = c(20, 40)),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 50L),
      metrics = list(rmse_test_metric()),
      n_replicates = 2L,
      seed = 1L
    )
    info <- preflight(config, condensed = FALSE)
    expect_equal(info$n_tasks, 4L)
    expect_equal(info$data_grid, 2L)
    expect_equal(info$fit_grid, 1L)
    expect_equal(info$n_replicates, 2L)
  })

  it("condensed mode runs without error", {
    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 50L),
      metrics = list(),
      n_replicates = 1L,
      seed = 1L
    )
    expect_invisible(preflight(config, condensed = TRUE))
  })
})

describe("F2 failed_tasks", {
  it("returns an empty tibble when all tasks succeed", {
    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 50L),
      metrics = list(),
      n_replicates = 2L,
      seed = 1L
    )
    res <- run_simulation(config, resume = "never", progress = FALSE)
    ft <- failed_tasks(res)
    expect_s3_class(ft, "data.frame")
    expect_equal(nrow(ft), 0L)
  })
})

describe("F3 as_tibble", {
  it("returns the summary tibble", {
    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 50L),
      metrics = list(posterior_summary_metric()),
      n_replicates = 2L,
      seed = 1L
    )
    res <- run_simulation(config, resume = "never", progress = FALSE)
    tb <- tibble::as_tibble(res)
    expect_s3_class(tb, "tbl_df")
    expect_equal(nrow(tb), 2L)
  })
})

describe("F4 seed error message", {
  it("errors with a helpful message when seed is missing", {
    expect_error(
      simulation_config(
        data_grid = data.frame(n = 20),
        fit_grid = data.frame(model = "lm"),
        data_generator = .gen,
        fitter = LinearRegressionFitter(n_draws = 50L)
      ),
      "RNG stream"
    )
  })
})
