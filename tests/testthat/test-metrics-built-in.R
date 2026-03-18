# test-metrics-built-in.R
# Tests for built-in metric classes: RmseMetric, BiasMetric, CoverageMetric, PosteriorMeanMetric

library(bayesim)

# Helper: create a valid fit result with draws
make_fit_result_with_draws <- function(draws = NULL, params = c(alpha = 1.0, beta = 2.0)) {
  if (is.null(draws)) {
    draws <- matrix(rnorm(500), ncol = 2, nrow = 250)
    colnames(draws) <- names(params)
  }
  new_fit_result(
    success = TRUE,
    draws = draws,
    diagnostics = list(rhat = 1.01),
    timing = list(total = 1.0)
  )
}

# Helper: create a valid data bundle for prediction metrics
make_data_bundle_for_predictions <- function(n = 50, true_alpha = 1.0, true_beta = 2.0) {
  set.seed(42)
  x <- rnorm(n)
  y <- true_alpha + true_beta * x + rnorm(n, sd = 0.5)
  train_idx <- seq_len(floor(n * 0.8))
  test_idx <- setdiff(seq_len(n), train_idx)
  list(
    train = data.frame(x = x[train_idx], y = y[train_idx]),
    test = data.frame(x = x[test_idx], y = y[test_idx]),
    response = "y",
    true_params = c(alpha = true_alpha, beta = true_beta),
    vars_of_interest = c("alpha", "beta"),
    references = c(alpha = 0, beta = 0),
    meta = list()
  )
}

# Helper: create context with predictions
make_prediction_context <- function(n_test, true_alpha = 1.0, true_beta = 2.0) {
  set.seed(123)
  # Predictions close to truth
  list(
    predictions = list(
      predicted_mean = true_alpha + true_beta * rnorm(n_test, mean = 1, sd = 0.3),
      predicted_samples = matrix(rnorm(n_test * 100), ncol = 100),
      predicted_sd = rep(0.5, n_test)
    )
  )
}

describe("RmseMetric", {
  it("creates valid metric with default name", {
    m <- rmse_metric()
    expect_true(inherits(m, "RmseMetric") || grepl("RmseMetric", class(m)[1]))
    expect_equal(m@name, "rmse")
    expect_equal(m@needs, "predictions")
    expect_false(m@required)
  })

  it("creates metric with custom name", {
    m <- rmse_metric(name = "my_rmse")
    expect_equal(m@name, "my_rmse")
  })

  it("computes RMSE correctly with test data", {
    m <- rmse_metric()
    db <- make_data_bundle_for_predictions(n = 100)
    n_test <- nrow(db$test)
    ctx <- make_prediction_context(n_test)

    output <- compute(m, make_fit_result_with_draws(), db, ctx, list(task_id = "test"))

    expect_true("value" %in% names(output))
    expect_true("n_obs" %in% names(output))
    expect_true(is.numeric(output$value))
    expect_true(output$value >= 0)
    expect_equal(output$n_obs, n_test)
  })

  it("returns NA when predictions are NULL", {
    m <- rmse_metric()
    db <- make_data_bundle_for_predictions()
    ctx <- list(predictions = NULL)

    output <- compute(m, make_fit_result_with_draws(), db, ctx, list(task_id = "test"))

    expect_equal(output$value, NA_real_)
  })

  it("returns NA when predictions context is missing", {
    m <- rmse_metric()
    db <- make_data_bundle_for_predictions()

    output <- compute(m, make_fit_result_with_draws(), db, list(), list(task_id = "test"))

    expect_equal(output$value, NA_real_)
  })

  it("falls back to train when test is NULL", {
    m <- rmse_metric()
    db <- make_data_bundle_for_predictions(n = 100)
    db$test <- NULL
    n_train <- nrow(db$train)
    ctx <- make_prediction_context(n_train)

    output <- compute(m, make_fit_result_with_draws(), db, ctx, list(task_id = "test"))

    expect_true("value" %in% names(output))
    expect_equal(output$n_obs, n_train)
  })
})

describe("BiasMetric", {
  it("creates valid metric with default name", {
    m <- bias_metric()
    expect_true(inherits(m, "BiasMetric") || grepl("BiasMetric", class(m)[1]))
    expect_equal(m@name, "bias")
  })

  it("computes bias correctly", {
    m <- bias_metric()
    db <- make_data_bundle_for_predictions(n = 100)
    n_test <- nrow(db$test)
    ctx <- make_prediction_context(n_test)

    output <- compute(m, make_fit_result_with_draws(), db, ctx, list(task_id = "test"))

    expect_true("value" %in% names(output))
    expect_true(is.numeric(output$value))
  })

  it("returns NA when predictions are NULL", {
    m <- bias_metric()
    db <- make_data_bundle_for_predictions()
    ctx <- list(predictions = NULL)

    output <- compute(m, make_fit_result_with_draws(), db, ctx, list(task_id = "test"))

    expect_equal(output$value, NA_real_)
  })
})

describe("CoverageMetric", {
  it("creates valid metric with default parameters", {
    m <- coverage_metric()
    expect_true(inherits(m, "CoverageMetric") || grepl("CoverageMetric", class(m)[1]))
    expect_equal(m@name, "coverage")
    expect_equal(m@needs, character())
    expect_equal(m@prob, 0.95)
  })

  it("creates metric with custom prob", {
    m <- coverage_metric(prob = 0.80)
    expect_equal(m@prob, 0.80)
  })

  it("computes coverage with matching draws and true params", {
    m <- coverage_metric()
    # Create draws centered at true params with small variance for high coverage
    set.seed(42)
    draws <- matrix(c(
      rnorm(500, mean = 1.0, sd = 0.1),
      rnorm(500, mean = 2.0, sd = 0.1)
    ), ncol = 2, nrow = 250)
    colnames(draws) <- c("alpha", "beta")

    fit_result <- make_fit_result_with_draws(draws = draws)
    db <- list(
      true_params = c(alpha = 1.0, beta = 2.0),
      vars_of_interest = c("alpha", "beta")
    )

    output <- compute(m, fit_result, db, list(), list(task_id = "test"))

    expect_true("mean" %in% names(output))
    expect_true("by_param" %in% names(output))
    expect_true(output$mean >= 0 && output$mean <= 1)
  })

  it("returns NA when draws are NULL", {
    m <- coverage_metric()
    fit_result <- new_fit_result(
      success = TRUE,
      draws = NULL,
      diagnostics = list(),
      timing = list(total = 1.0)
    )
    db <- list(true_params = c(alpha = 1.0), vars_of_interest = "alpha")

    output <- compute(m, fit_result, db, list(), list(task_id = "test"))

    expect_equal(output$value, NA_real_)
  })

  it("returns NA when true_params are NULL", {
    m <- coverage_metric()
    fit_result <- make_fit_result_with_draws()
    db <- list(true_params = NULL, vars_of_interest = "alpha")

    output <- compute(m, fit_result, db, list(), list(task_id = "test"))

    expect_equal(output$value, NA_real_)
  })

  it("handles missing parameter in draws gracefully", {
    m <- coverage_metric()
    draws <- matrix(rnorm(500), ncol = 1, nrow = 500)
    colnames(draws) <- "alpha"

    fit_result <- make_fit_result_with_draws(draws = draws)
    db <- list(
      true_params = c(alpha = 1.0, beta = 2.0),
      vars_of_interest = c("alpha", "beta")
    )

    output <- compute(m, fit_result, db, list(), list(task_id = "test"))

    # beta should be NA since it's not in draws
    expect_true(is.na(output$by_param["beta"]))
    expect_true(!is.na(output$by_param["alpha"]))
  })
})

describe("PosteriorMeanMetric", {
  it("creates valid metric with default name", {
    m <- posterior_mean_metric()
    expect_true(inherits(m, "PosteriorMeanMetric") || grepl("PosteriorMeanMetric", class(m)[1]))
    expect_equal(m@name, "posterior_mean")
    expect_equal(m@needs, character())
  })

  it("computes posterior means correctly", {
    m <- posterior_mean_metric()
    # Use draws with known means
    set.seed(42)
    draws <- matrix(c(
      rep(1.0, 500),
      rep(2.0, 500)
    ), ncol = 2, nrow = 500)
    colnames(draws) <- c("alpha", "beta")

    fit_result <- make_fit_result_with_draws(draws = draws)
    db <- list(vars_of_interest = c("alpha", "beta"))

    output <- compute(m, fit_result, db, list(), list(task_id = "test"))

    expect_true("alpha" %in% names(output))
    expect_true("beta" %in% names(output))
    # With constant draws, posterior means should be exact
    expect_equal(output$alpha, 1.0)
    expect_equal(output$beta, 2.0)
  })

  it("returns NA when draws are NULL", {
    m <- posterior_mean_metric()
    fit_result <- new_fit_result(
      success = TRUE,
      draws = NULL,
      diagnostics = list(),
      timing = list(total = 1.0)
    )

    output <- compute(m, fit_result, list(vars_of_interest = "alpha"), list(), list(task_id = "test"))

    expect_equal(output$value, NA_real_)
  })

  it("uses colnames when vars_of_interest is NULL", {
    m <- posterior_mean_metric()
    set.seed(42)
    draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
    colnames(draws) <- c("a", "b")

    fit_result <- make_fit_result_with_draws(draws = draws)

    output <- compute(m, fit_result, list(), list(), list(task_id = "test"))

    expect_true("a" %in% names(output))
    expect_true("b" %in% names(output))
  })
})
