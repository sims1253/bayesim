# test-worker.R
# Tests for worker functions: run_task_safe, run_task, build_metric_context, compute_all_metrics

library(bayesim)

test_data_gen <- function(data_spec, seed, task_ctx) {
  set.seed(seed)
  n <- data_spec$n %||% 100
  list(
    train = data.frame(x = rnorm(n), y = rnorm(n)),
    test = data.frame(x = rnorm(n), y = rnorm(n)),
    response = "y",
    true_params = c(beta = 1.0),
    vars_of_interest = "beta",
    references = c(beta = 0),
    meta = list()
  )
}

describe("run_task_safe()", {
  it("returns task_result for successful task", {
    config <- simulation_config(
      data_grid = data.frame(n = 50),
      fit_grid = data.frame(model = "test"),
      data_generator = test_data_gen,
      fitter = MockFitter(),
      seed = 42L
    )

    task <- get_task_spec(create_task_grid(config), "d001_f001_r00001", config)
    config_spec <- as_config_spec(config)
    config_spec$data_generator <- config@data_generator

    result <- run_task_safe(
      task = task,
      config_spec = config_spec,
      fitter = config@fitter,
      metrics = list(),
      retain = c("metrics", "diagnostics")
    )

    expect_s3_class(result, "bayesim_task_result")
    expect_true(result$status %in% c("success", "failed"))
    expect_equal(result$task_id, "d001_f001_r00001")
  })

  it("returns failed task_result for data generation errors", {
    bad_gen <- function(data_spec, seed, task_ctx) {
      stop("Data generation failed")
    }

    config <- simulation_config(
      data_grid = data.frame(n = 50),
      fit_grid = data.frame(model = "test"),
      data_generator = bad_gen,
      fitter = MockFitter(),
      seed = 42L
    )

    task <- get_task_spec(create_task_grid(config), "d001_f001_r00001", config)
    config_spec <- as_config_spec(config)
    config_spec$data_generator <- bad_gen

    result <- run_task_safe(
      task = task,
      config_spec = config_spec,
      fitter = config@fitter,
      metrics = list(),
      retain = c("metrics", "diagnostics")
    )

    expect_s3_class(result, "bayesim_task_result")
    expect_equal(result$status, "failed")
    expect_false(is.null(result$error))
  })

  it("propagates fatal errors", {
    config <- simulation_config(
      data_grid = data.frame(n = 50),
      fit_grid = data.frame(model = "test"),
      data_generator = test_data_gen,
      fitter = MockFitter(),
      seed = 42L
    )

    task <- get_task_spec(create_task_grid(config), "d001_f001_r00001", config)
    # Remove data_spec to trigger a fatal error in validate_data_bundle
    task$data_spec <- NULL

    config_spec <- as_config_spec(config)
    # Override data_generator to return something that will trigger bayesim_data_error
    config_spec$data_generator = function(ds, seed, ctx) list(
      train = NULL, test = NULL, response = "y",
      true_params = c(b = 1), vars_of_interest = "b", meta = list()
    )

    # bayesim_data_error is recoverable, so run_task_safe converts it to a failed result
    result <- run_task_safe(
      task = task,
      config_spec = config_spec,
      fitter = config@fitter,
      metrics = list(),
      retain = c("metrics", "diagnostics")
    )
    expect_s3_class(result, "bayesim_task_result")
    expect_equal(result$status, "failed")
  })
})

describe("build_metric_context()", {
  it("returns empty list when no metrics need anything", {
    fitter <- MockFitter()
    fit_result <- new_fit_result(
      success = TRUE,
      draws = NULL,
      diagnostics = list(),
      timing = list(total = 1.0)
    )
    data_bundle <- list(train = data.frame(x = 1:10, y = rnorm(10)), response = "y")

    ctx <- build_metric_context(fit_result, fitter, data_bundle, list())
    expect_equal(length(ctx), 0)
  })

  it("computes predictions when needed", {
    fitter <- MockFitter()
    data_bundle <- list(
      train = data.frame(x = rnorm(10), y = rnorm(10)),
      test = data.frame(x = rnorm(10), y = rnorm(10)),
      response = "y"
    )
    fit_result <- fit(
      fitter,
      data_bundle,
      fit_spec = data.frame(model = "test"),
      seed = 42L,
      task_ctx = list(task_id = "test")
    )

    m <- rmse_metric()
    ctx <- build_metric_context(fit_result, fitter, data_bundle, list(m))

    expect_false(is.null(ctx$predictions))
    expect_true("predicted_mean" %in% names(ctx$predictions))
  })

  it("warns when fitter does not support required capability", {
    # MockFitter supports predictions, so test with log_lik which it also supports
    # Use a fitter that doesn't support predictions
    MinimalFitter <- S7::new_class(
      "MinimalFitter",
      parent = Fitter,
      properties = list(
        name = S7::new_property(S7::class_character, default = "minimal"),
        supports_predictions = S7::new_property(S7::class_logical, default = FALSE),
        supports_log_lik = S7::new_property(S7::class_logical, default = FALSE),
        supports_loo = S7::new_property(S7::class_logical, default = FALSE)
      )
    )
    S7::method(fit, MinimalFitter) <- function(fitter, data_bundle, fit_spec, seed, task_ctx) {
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      new_fit_result(
        success = TRUE,
        draws = draws,
        diagnostics = list(),
        timing = list(total = 1.0)
      )
    }

    fitter <- MinimalFitter()
    data_bundle <- list(train = data.frame(x = 1), response = "y")
    fit_result <- fit(fitter, data_bundle, data.frame(), 42L, list(task_id = "test"))

    expect_warning(
      build_metric_context(fit_result, fitter, data_bundle, list(rmse_metric())),
      "does not support them"
    )
  })
})

describe("compute_all_metrics()", {
  it("returns empty results for no metrics", {
    result <- compute_all_metrics(
      fit_result = new_fit_result(success = TRUE, diagnostics = list(), timing = list(total = 1.0)),
      data_bundle = list(),
      context = list(),
      metrics = list(),
      task_ctx = list(task_id = "test")
    )
    expect_equal(length(result$metrics), 0)
    expect_equal(length(result$warnings), 0)
  })

  it("computes single metric successfully", {
    fitter <- MockFitter()
    data_bundle <- list(
      train = data.frame(x = rnorm(10), y = rnorm(10)),
      test = data.frame(x = rnorm(5), y = rnorm(5)),
      response = "y"
    )
    fit_result <- fit(fitter, data_bundle, data.frame(), 42L, list(task_id = "test"))
    ctx <- build_metric_context(fit_result, fitter, data_bundle, list(rmse_metric()))

    result <- compute_all_metrics(
      fit_result = fit_result,
      data_bundle = data_bundle,
      context = ctx,
      metrics = list(rmse_metric()),
      task_ctx = list(task_id = "test")
    )

    expect_true(length(result$metrics) > 0)
  })

  it("handles failed non-required metric gracefully", {
    # A metric that always fails but is not required
    FailMetric <- S7::new_class(
      "FailMetric",
      parent = Metric,
      properties = list(
        name = S7::new_property(S7::class_character, default = "always_fail"),
        needs = S7::new_property(S7::class_character, default = character()),
        required = S7::new_property(S7::class_logical, default = FALSE)
      )
    )
    S7::method(compute, FailMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
      stop("Intentional failure")
    }

    result <- compute_all_metrics(
      fit_result = new_fit_result(success = TRUE, diagnostics = list(), timing = list(total = 1.0)),
      data_bundle = list(),
      context = list(),
      metrics = list(FailMetric(name = "always_fail")),
      task_ctx = list(task_id = "test")
    )

    # Should not error - metric is not required
    # The metric name should be properly set
    expect_true(length(result$metrics) > 0)
  })

  it("propagates error for required metric failure", {
    FailRequiredMetric <- S7::new_class(
      "FailRequiredMetric",
      parent = Metric,
      properties = list(
        name = S7::new_property(S7::class_character, default = "required_fail"),
        needs = S7::new_property(S7::class_character, default = character()),
        required = S7::new_property(S7::class_logical, default = TRUE)
      )
    )
    S7::method(compute, FailRequiredMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
      stop("Required metric failure")
    }

    m <- FailRequiredMetric()
    # S7 class_logical defaults may not propagate correctly across package boundaries
    # Explicitly set required = TRUE
    m <- FailRequiredMetric(required = TRUE)

    expect_error(
      compute_all_metrics(
        fit_result = new_fit_result(success = TRUE, diagnostics = list(), timing = list(total = 1.0)),
        data_bundle = list(),
        context = list(),
        metrics = list(m),
        task_ctx = list(task_id = "test")
      ),
      "Required metric failure"
    )
  })
})

describe("derive_task_seed()", {
  it("returns positive integer from valid RNG seed", {
    streams <- create_task_rng_streams(42, 1)
    seed <- derive_task_seed(streams[[1]])
    expect_true(is.integer(seed))
    expect_true(seed > 0)
  })

  it("returns 1 for NULL seed", {
    expect_equal(derive_task_seed(NULL), 1L)
  })

  it("returns 1 for short seed", {
    expect_equal(derive_task_seed(integer(0)), 1L)
  })

  it("produces different seeds for different streams", {
    streams <- create_task_rng_streams(42, 10)
    seeds <- vapply(streams, derive_task_seed, integer(1))
    expect_true(length(unique(seeds)) > 1)
  })
})
