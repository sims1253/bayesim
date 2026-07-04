# tests/testthat/test-parquet-summary.R
# Workstream I8: optional parquet summary sidecar.

skip_on_cran()

library(bayesim)
library(testthat)

# Tiny data generator producing a valid data_bundle (y = beta*x + noise).
# Mirrors bayesim_example_data_generator() shape so LinearRegressionFitter works.
.datagen <- function(data_spec, task_ctx) {
  n <- as.integer(data_spec$n %||% 20L)
  beta <- as.numeric(data_spec$beta %||% 1)
  sigma <- as.numeric(data_spec$sigma %||% 1)
  x <- stats::rnorm(n)
  y <- beta * x + stats::rnorm(n, sd = sigma)
  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(beta = beta, sigma = sigma),
    vars_of_interest = c("beta", "sigma")
  )
}

.fitter <- bayesim::LinearRegressionFitter(n_draws = 50L)
.metric <- bayesim::pred_rmse_metric()

describe("summary_format = 'parquet'", {
  it("writes summary.parquet next to the checkpoint and round-trips", {
    skip_if_not(
      requireNamespace("nanoparquet", quietly = TRUE),
      message = "nanoparquet not available"
    )

    tmp <- tempfile()
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

    data_grid <- data.frame(n = 20L, beta = 1, sigma = 1)
    fit_grid <- data.frame(model = "baseline", stringsAsFactors = FALSE)

    config <- simulation_config(
      data_grid = data_grid,
      fit_grid = fit_grid,
      data_generator = .datagen,
      fitter = .fitter,
      metrics = list(.metric),
      n_replicates = 2L,
      seed = 1L,
      result_path = tmp,
      summary_format = "parquet"
    )

    expect_true(bayesim:::is_simulation_config(config))
    expect_equal(config@summary_format, "parquet")

    result <- run_simulation(config, progress = FALSE)

    parquet_path <- file.path(tmp, "summary.parquet")
    expect_true(file.exists(parquet_path))

    # Round-trip via read_summary().
    df <- read_summary(parquet_path)
    expect_s3_class(df, "data.frame")
    expect_gt(nrow(df), 0L)
    # Core columns are present.
    expect_true("task_id" %in% names(df))
    expect_true("status" %in% names(df))

    # read_summary dispatches by extension on an rds file.
    rds_path <- file.path(tempdir(), "bayesim_summary_test.rds")
    saveRDS(data.frame(a = 1:3, b = letters[1:3]), rds_path)
    on.exit(unlink(rds_path), add = TRUE)
    expect_equal(nrow(read_summary(rds_path)), 3L)
  })

  it("defaults to 'rds' and writes no summary.parquet", {
    skip_if_not(
      requireNamespace("nanoparquet", quietly = TRUE),
      message = "nanoparquet not available"
    )

    tmp <- tempfile()
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

    config <- simulation_config(
      data_grid = data.frame(n = 20L, beta = 1, sigma = 1),
      fit_grid = data.frame(model = "baseline", stringsAsFactors = FALSE),
      data_generator = .datagen,
      fitter = .fitter,
      metrics = list(.metric),
      n_replicates = 1L,
      seed = 2L,
      result_path = tmp
    )

    expect_equal(config@summary_format, "rds")
    run_simulation(config, progress = FALSE)
    expect_false(file.exists(file.path(tmp, "summary.parquet")))
  })

  it("rejects invalid summary_format values", {
    skip_if_not(
      requireNamespace("nanoparquet", quietly = TRUE),
      message = "nanoparquet not available"
    )
    expect_error(
      simulation_config(
        data_grid = data.frame(n = 20L),
        fit_grid = data.frame(model = "baseline"),
        data_generator = .datagen,
        fitter = .fitter,
        metrics = list(.metric),
        seed = 3L,
        summary_format = "csv"
      )
    )
  })
})
