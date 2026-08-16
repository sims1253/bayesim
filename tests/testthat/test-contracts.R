# test-contracts.R
# Contract tests for the core infrastructure: error classes, result
# constructors, the fitter/metric seams, config canonicalization, and schema
# validation.
#
# Per the testing strategy in IMPROVEMENT_PLAN.md this file intentionally does
# NOT re-assert S7 class hierarchies, constructor defaults, or every
# single-argument invalid permutation already covered by property validators.
# It keeps the contracts that behavior elsewhere depends on: error
# classification (fatal vs recoverable), task-result invariants, metric
# output schema validation, fingerprint canonicalization, and data-bundle
# schema validation.

# Load the package
library(bayesim)

# =============================================================================
# 1. Error Classes (from errors.R)
# =============================================================================

describe("Error Classes", {
  it("constructs every error type with its class, message, and base classes", {
    errors <- list(
      bayesim_error = bayesim_error,
      bayesim_config_error = bayesim_config_error,
      bayesim_contract_error = bayesim_contract_error,
      bayesim_checkpoint_error = bayesim_checkpoint_error,
      bayesim_internal_error = bayesim_internal_error,
      bayesim_data_error = bayesim_data_error,
      bayesim_fit_error = bayesim_fit_error,
      bayesim_metric_error = bayesim_metric_error
    )
    for (error_name in names(errors)) {
      err <- errors[[error_name]]("Test error message")
      expect_s3_class(err, error_name)
      expect_s3_class(err, "bayesim_error")
      expect_s3_class(err, "error")
      expect_equal(err$message, "Test error message")
      expect_true(is_bayesim_error(err))
    }
    expect_false(is_bayesim_error(simpleError("Test")))
    expect_false(is_bayesim_error(list()))
  })

  it("partitions fatal vs recoverable errors", {
    # The engine's error handling depends on this exact partition: fatal
    # conditions stop the run, recoverable ones become failed task results.
    fatal <- list(
      bayesim_config_error,
      bayesim_contract_error,
      bayesim_checkpoint_error,
      bayesim_internal_error
    )
    recoverable <- list(
      bayesim_data_error,
      bayesim_fit_error,
      bayesim_metric_error
    )
    for (ctor in fatal) {
      expect_true(is_fatal_error(ctor("Test")), info = class(ctor)[1])
      expect_false(is_recoverable_error(ctor("Test")), info = class(ctor)[1])
    }
    for (ctor in recoverable) {
      expect_false(is_fatal_error(ctor("Test")), info = class(ctor)[1])
      expect_true(is_recoverable_error(ctor("Test")), info = class(ctor)[1])
    }
    # The base class and foreign errors are neither.
    expect_false(is_fatal_error(bayesim_error("Test")))
    expect_false(is_recoverable_error(bayesim_error("Test")))
    expect_false(is_fatal_error(simpleError("Test")))
    expect_false(is_recoverable_error(simpleError("Test")))
  })

  it("dispatches through the bayesim_error class in tryCatch handlers", {
    result <- tryCatch(
      stop(bayesim_config_error("Config error")),
      bayesim_error = function(e) {
        if (is_fatal_error(e)) "fatal" else "recoverable"
      }
    )
    expect_equal(result, "fatal")
  })
})

# =============================================================================
# 2. Result Constructors (from results.R)
# =============================================================================

describe("Result Constructors", {
  describe("new_fit_result()", {
    it("round-trips a successful fit with all fields", {
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")

      result <- new_fit_result(
        success = TRUE,
        fit = list(mock = TRUE),
        draws = draws,
        diagnostics = list(rhat = c(alpha = 1.01, beta = 1.00)),
        timing = list(total = 10.5, warmup = 5.0, sample = 5.5),
        warnings = c("Warning 1", "Warning 2")
      )

      expect_s3_class(result, "bayesim_fit_result")
      expect_true(result$success)
      expect_equal(dim(result$draws), c(50, 2))
      expect_equal(colnames(result$draws), c("alpha", "beta"))
      expect_equal(length(result$warnings), 2)
      expect_equal(result$timing$total, 10.5)
    })

    it("round-trips a failed fit with its error condition", {
      err <- bayesim_fit_error("Convergence failed")
      result <- new_fit_result(
        success = FALSE,
        error = err,
        timing = list(total = 2.0)
      )

      expect_false(result$success)
      expect_identical(result$error, err)
    })

    it("enforces the success/error consistency invariants", {
      expect_error(
        new_fit_result(
          success = TRUE,
          error = simpleError("Error"),
          timing = list(total = 1.0)
        ),
        "When success is TRUE, error must be NULL"
      )
      expect_error(
        new_fit_result(
          success = FALSE,
          error = NULL,
          timing = list(total = 1.0)
        ),
        "When success is FALSE, error must be non-NULL"
      )
    })

    it("requires draws to be a named matrix (the metric-pipeline contract)", {
      expect_error(
        new_fit_result(
          success = TRUE,
          draws = matrix(rnorm(100), ncol = 2),
          timing = list(total = 1.0)
        ),
        "draws matrix must have colnames"
      )
      expect_error(
        new_fit_result(
          success = TRUE,
          draws = data.frame(a = 1:10),
          timing = list(total = 1.0)
        ),
        "draws must be a matrix when not NULL"
      )
    })
  })

  describe("new_task_result()", {
    it("creates valid success, failed, and skipped task results", {
      success <- new_task_result(
        task_id = "d001_f001_r00001",
        status = "success",
        metrics = list(rmse = 0.05),
        diagnostics = list(n_eff = 500),
        timing = list(total = 5.2)
      )
      expect_s3_class(success, "bayesim_task_result")
      expect_equal(success$task_id, "d001_f001_r00001")
      expect_equal(success$status, "success")

      failed <- new_task_result(
        task_id = "d001_f001_r00002",
        status = "failed",
        error = list(
          error_class = "convergence_error",
          error_message = "R-hat > 1.1"
        ),
        timing = list(total = 2.0)
      )
      expect_equal(failed$status, "failed")
      expect_equal(failed$error$error_class, "convergence_error")

      skipped <- new_task_result(
        task_id = "d001_f001_r00003",
        status = "skipped",
        timing = list(total = 0.0)
      )
      expect_equal(skipped$status, "skipped")
    })

    it("enforces the status/payload invariants", {
      expect_error(
        new_task_result(
          task_id = "test",
          status = "success",
          metrics = NULL,
          timing = list(total = 1.0)
        ),
        "When status is 'success', metrics must not be NULL"
      )
      expect_error(
        new_task_result(
          task_id = "test",
          status = "failed",
          error = NULL,
          timing = list(total = 1.0)
        ),
        "When status is 'failed', error must not be NULL"
      )
    })
  })

  describe("new_simulation_result()", {
    it("creates a valid simulation result from task results", {
      task1 <- new_task_result(
        task_id = "d001_f001_r00001",
        status = "success",
        metrics = list(rmse = 0.05),
        timing = list(total = 5.0)
      )

      result <- new_simulation_result(
        config_fingerprint = "abc123",
        task_results = list(task1),
        summary = data.frame(task_id = "d001_f001_r00001", rmse = 0.05),
        timing = list(total = 10.0),
        errors = data.frame(task_id = character(), error_message = character())
      )

      expect_s3_class(result, "bayesim_simulation_result")
      expect_equal(result$config_fingerprint, "abc123")
      expect_equal(length(result$task_results), 1L)
    })

    it("requires every element of task_results to be a task result", {
      expect_error(
        new_simulation_result(
          config_fingerprint = "test",
          task_results = list(list(not_a_result = TRUE))
        ),
        "All elements of task_results must be bayesim_task_result objects"
      )
    })

    it("materializes NULL inputs into usable empty defaults", {
      result <- new_simulation_result(
        config_fingerprint = "test",
        task_results = NULL,
        summary = NULL,
        timing = NULL,
        errors = NULL
      )
      expect_equal(result$task_results, list())
      expect_true(is.data.frame(result$summary))
      expect_equal(result$timing$total, 0)
      expect_true(is.data.frame(result$errors))
    })
  })
})

# =============================================================================
# 3. Fitter Class (from fitter.R)
# =============================================================================

describe("Fitter Class", {
  # MockFitter is the reference implementation of the fitter seam; these
  # checks define the shapes every fitter must produce.
  describe("MockFitter seam contracts", {
    .fit_once <- function(fitter = MockFitter()) {
      data_bundle <- list(
        train = data.frame(x = 1:10, y = rnorm(10)),
        test = NULL,
        response = "y"
      )
      list(
        bundle = data_bundle,
        fit = fit_model(
          fitter,
          data_bundle,
          fit_spec = data.frame(model = "test"),
          seed = 42L,
          task_ctx = list(task_id = "test")
        )
      )
    }

    it("fit_model() returns a bayesim_fit_result", {
      fitter <- MockFitter()
      out <- .fit_once(fitter)
      expect_s3_class(out$fit, "bayesim_fit_result")
      expect_equal(out$fit$fit$fitter, "mock")
    })

    it("extract_draws() returns a named matrix, respecting variables", {
      fitter <- MockFitter()
      out <- .fit_once(fitter)

      draws <- extract_draws(fitter, out$fit)
      expect_true(is.matrix(draws))
      expect_true(all(c("intercept", "beta", "sigma") %in% colnames(draws)))

      sel <- extract_draws(fitter, out$fit, variables = c("beta", "sigma"))
      expect_equal(colnames(sel), c("beta", "sigma"))
    })

    it("predict_fit() returns mean/sd/samples in S x N orientation", {
      fitter <- MockFitter()
      out <- .fit_once(fitter)

      preds <- predict_fit(fitter, out$fit)
      expect_setequal(
        names(preds),
        c("predicted_mean", "predicted_samples", "predicted_sd")
      )
      expect_equal(length(preds$predicted_mean), 10)
      expect_equal(ncol(preds$predicted_samples), 10)
    })

    it("log_lik_matrix() returns an S x N matrix", {
      fitter <- MockFitter()
      out <- .fit_once(fitter)

      ll <- log_lik_matrix(fitter, out$fit)
      expect_true(is.matrix(ll))
      expect_equal(ncol(ll), 10)
      expect_equal(
        nrow(ll),
        as.integer(fitter@n_draws) * as.integer(fitter@n_chains)
      )
    })

    it("loo_fit() and fit_diagnostics() return documented structures", {
      fitter <- MockFitter()
      out <- .fit_once(fitter)

      loo_result <- loo_fit(fitter, out$fit)
      expect_setequal(
        names(loo_result),
        c("elpd", "p_loo", "elpd_se", "pareto_k")
      )

      diag <- fit_diagnostics(fitter, out$fit)
      expect_true(all(c("rhat_max", "ess_bulk", "divergent") %in% names(diag)))
    })
  })

  it("supports custom Fitter subclasses via S7 method registration", {
    TestFitter <- S7::new_class(
      "TestFitter",
      parent = Fitter,
      properties = list(
        name = S7::new_property(S7::class_character, default = "test")
      )
    )
    S7::method(fit_model, TestFitter) <- function(
      fitter,
      data_bundle,
      fit_spec,
      seed,
      task_ctx
    ) {
      draws <- matrix(rnorm(40), ncol = 2)
      colnames(draws) <- c("a", "b")
      new_fit_result(
        success = TRUE,
        draws = draws,
        timing = list(total = 1.0)
      )
    }

    fitter <- TestFitter()
    expect_equal(fitter@name, "test")

    fit_result <- fit_model(
      fitter,
      data.frame(x = 1),
      data.frame(),
      1L,
      list(task_id = "test")
    )
    expect_s3_class(fit_result, "bayesim_fit_result")

    # Unimplemented seam methods fail loudly rather than silently degrading.
    expect_error(
      extract_draws(fitter, fit_result),
      "Can't find method"
    )
  })
})

# =============================================================================
# 4. Metric Class (from metric.R)
# =============================================================================

describe("Metric Class", {
  describe("validate_metric_output()", {
    it("accepts valid scalar, vector, and mixed outputs", {
      expect_silent(validate_metric_output(
        list(rmse = 0.5, n_obs = 100L),
        "test_metric"
      ))
      expect_silent(validate_metric_output(
        list(params = c(alpha = 0.1, beta = 0.2, gamma = 0.3)),
        "test_metric"
      ))
      expect_silent(validate_metric_output(
        list(rmse = 0.5, params = c(alpha = 0.1, beta = 0.2)),
        "test_metric"
      ))
    })

    it("rejects structurally invalid outputs", {
      expect_error(
        validate_metric_output(c(a = 1, b = 2), "test_metric"),
        "must be a list"
      )
      expect_error(
        validate_metric_output(list(), "test_metric"),
        "cannot be empty"
      )
      expect_error(
        validate_metric_output(list(0.5, 1.0), "test_metric"),
        "unnamed or empty-named"
      )
      expect_error(
        validate_metric_output(list(rmse = NULL), "test_metric"),
        "is NULL \\(not allowed\\)"
      )
      expect_error(
        validate_metric_output(list(nested = list(a = 1)), "test_metric"),
        "must be either"
      )
      expect_error(
        validate_metric_output(list(df = data.frame(a = 1:3)), "test_metric"),
        "must be either"
      )
      expect_error(
        validate_metric_output(list(mat = matrix(1:4, 2, 2)), "test_metric"),
        "must be either"
      )
      expect_error(
        validate_metric_output(list(params = c(1, 2, 3)), "test_metric"),
        "unnamed numeric vector"
      )
    })
  })

  describe("flatten_metric_output()", {
    it("prefixes scalars and named vectors with the metric name", {
      flat <- flatten_metric_output(
        list(rmse = 0.5, params = c(alpha = 0.1, beta = 0.2)),
        "my_metric"
      )

      expect_equal(flat$my_metric__rmse, 0.5)
      expect_equal(flat$my_metric__params__alpha, 0.1)
      expect_equal(flat$my_metric__params__beta, 0.2)
    })

    it("validates input before flattening", {
      expect_error(
        flatten_metric_output(list(unnamed = NULL), "test"),
        "is NULL"
      )
    })
  })

  it("supports custom Metric subclasses via S7 method registration", {
    RMSEMetric <- S7::new_class(
      "RMSEMetric",
      parent = Metric,
      properties = list(
        name = S7::new_property(S7::class_character, default = "rmse"),
        needs = S7::new_property(
          S7::class_character,
          default = "predictions"
        ),
        required = S7::new_property(S7::class_logical, default = FALSE)
      )
    )
    S7::method(compute_metric, RMSEMetric) <- function(
      metric,
      fit_result,
      data_bundle,
      context,
      task_ctx
    ) {
      preds <- context$predictions$predicted_mean
      actual <- data_bundle$test[[data_bundle$response]]
      list(
        value = sqrt(mean((preds - actual)^2)),
        n_obs = length(actual)
      )
    }

    metric <- RMSEMetric()
    output <- compute_metric(
      metric,
      list(),
      list(test = data.frame(y = c(1, 2, 3)), response = "y"),
      list(predictions = list(predicted_mean = c(1.1, 1.9, 3.2))),
      list(task_id = "test")
    )

    expect_silent(validate_metric_output(output, "rmse"))
    expect_setequal(names(output), c("value", "n_obs"))

    # Unimplemented methods fail loudly.
    BareMetric <- S7::new_class(
      "BareMetric",
      parent = Metric,
      properties = list(
        name = S7::new_property(S7::class_character, default = "bare")
      )
    )
    expect_error(
      compute_metric(BareMetric(), NULL, NULL, NULL, NULL),
      "Can't find method"
    )
  })
})

# =============================================================================
# 5. SimulationConfig (from simulation-config.R)
# =============================================================================

describe("SimulationConfig", {
  test_data_gen <- function(data_spec, task_ctx) {
    list(
      train = data.frame(x = 1:10, y = rnorm(10)),
      test = NULL,
      response = "y"
    )
  }

  describe("simulation_config()", {
    it("creates valid objects with the full parameter surface", {
      config <- simulation_config(
        data_grid = data.frame(n = c(100, 500)),
        fit_grid = data.frame(model = c("A", "B")),
        data_generator = test_data_gen,
        fitter = MockFitter(),
        metrics = list(),
        n_replicates = 100L,
        seed = 42L,
        result_path = "/tmp/results",
        checkpoint_format = "rds",
        checkpoint_every = 25L,
        retain = list(success = c("metrics", "diagnostics"), error = "debug"),
        max_errors = 10
      )

      expect_true(is_simulation_config(config))
      expect_equal(config@n_replicates, 100L)
      expect_equal(config@seed, 42L)
      expect_equal(config@result_path, "/tmp/results")
      expect_equal(config@checkpoint_every, 25L)
      expect_equal(config@retain$success, c("metrics", "diagnostics"))
      expect_true(all(c("metrics", "draws", "fit") %in% config@retain$error))
    })

    it("accepts explicit task_grid specifications", {
      config <- simulation_config(
        task_grid = tibble::tibble(
          data_spec = list(list(n = 50), list(n = 100)),
          fit_spec = list(list(model = "a"), list(model = "b")),
          rep_idx = c(1L, 2L)
        ),
        data_generator = test_data_gen,
        seed = 42L
      )

      expect_true(is_simulation_config(config))
      expect_equal(nrow(config@task_grid), 2)
      expect_null(config@data_grid)
      expect_null(config@fit_grid)
    })

    it("validates the retain option set and metric object types", {
      expect_error(
        simulation_config(
          data_grid = data.frame(a = 1),
          fit_grid = data.frame(a = 1),
          data_generator = test_data_gen,
          retain = c("invalid_option"),
          seed = 42L
        ),
        "retain contains invalid options"
      )
      expect_error(
        simulation_config(
          data_grid = data.frame(a = 1),
          fit_grid = data.frame(a = 1),
          data_generator = test_data_gen,
          metrics = c("rmse", "bias"),
          seed = 42L
        ),
        "metrics must be Metric objects"
      )
    })
  })

  describe("as_config_spec()", {
    it("produces the design spec and excludes runtime policy", {
      config <- simulation_config(
        data_grid = data.frame(n = 100),
        fit_grid = data.frame(model = "baseline"),
        data_generator = test_data_gen,
        result_path = "/tmp/test",
        checkpoint_every = 100L,
        seed = 42L
      )

      spec <- as_config_spec(config)

      # Design fields enter the canonical spec...
      expect_true(all(
        c(
          "data_grid",
          "fit_grid",
          "data_generator_spec",
          "seed"
        ) %in%
          names(spec)
      ))
      # ...runtime-specific fields are excluded (fingerprint canonicalization).
      expect_false("result_path" %in% names(spec))
      expect_false("checkpoint_every" %in% names(spec))
      expect_false("keep_checkpoints" %in% names(spec))

      expect_error(as_config_spec(list()), "config must be a SimulationConfig")
    })
  })

  describe("compute_config_fingerprint()", {
    it("is stable, discriminating, and SHA-256 shaped", {
      config <- simulation_config(
        data_grid = data.frame(n = 100),
        fit_grid = data.frame(model = "baseline"),
        data_generator = test_data_gen,
        seed = 42L
      )
      variant <- simulation_config(
        data_grid = data.frame(n = 200),
        fit_grid = data.frame(model = "baseline"),
        data_generator = test_data_gen,
        seed = 42L
      )

      fp <- compute_config_fingerprint(config)
      expect_equal(fp, compute_config_fingerprint(config))
      expect_equal(nchar(fp), 64L)
      expect_false(isTRUE(all.equal(fp, compute_config_fingerprint(variant))))

      expect_error(
        compute_config_fingerprint(list()),
        "config must be a SimulationConfig object"
      )
    })

    it("B4: excludes runtime-only fields from study identity", {
      base <- list(
        data_grid = data.frame(a = 1),
        fit_grid = data.frame(a = 1),
        data_generator = test_data_gen,
        seed = 42L,
        n_replicates = 5L
      )
      cfg1 <- do.call(simulation_config, c(base, list(retain = c("metrics"))))
      cfg2 <- do.call(
        simulation_config,
        c(base, list(retain = c("metrics", "diagnostics")))
      )
      expect_equal(
        compute_config_fingerprint(cfg1),
        compute_config_fingerprint(cfg2)
      )
      expect_identical(
        config_fingerprint(cfg1),
        compute_config_fingerprint(cfg1)
      )
    })
  })

  describe("get_total_tasks()", {
    it("multiplies grid sizes by replicates", {
      config <- simulation_config(
        data_grid = data.frame(n = c(100, 500)),
        fit_grid = data.frame(model = c("A", "B")),
        data_generator = test_data_gen,
        n_replicates = 100L,
        seed = 42L
      )
      expect_equal(get_total_tasks(config), 2 * 2 * 100)

      single <- simulation_config(
        data_grid = data.frame(n = 100),
        fit_grid = data.frame(model = "A"),
        data_generator = test_data_gen,
        n_replicates = 50L,
        seed = 42L
      )
      expect_equal(get_total_tasks(single), 50)
    })

    it("counts explicit task_grid rows directly", {
      config <- simulation_config(
        task_grid = tibble::tibble(
          data_spec = list(list(n = 10), list(n = 20), list(n = 30)),
          fit_spec = list(
            list(model = "a"),
            list(model = "a"),
            list(model = "b")
          )
        ),
        data_generator = test_data_gen,
        seed = 42L
      )
      expect_equal(get_total_tasks(config), 3)
      expect_error(get_total_tasks(list()), "config must be a SimulationConfig")
    })
  })
})

# =============================================================================
# 6. Utility Functions (from utils.R)
# =============================================================================

describe("Utility Functions", {
  describe("capture_error_info()", {
    it("captures class, message, call, and traceback from conditions", {
      err <- tryCatch(stop("Test error message"), error = function(e) e)
      info <- capture_error_info(err)

      expect_true(all(
        c(
          "error_class",
          "error_message",
          "call",
          "traceback"
        ) %in%
          names(info)
      ))
      expect_equal(info$error_message, "Test error message")
      expect_lte(length(info$traceback), 20)
    })

    it("preserves calls from the original error site", {
      inner_failure <- function() stop("nested error")
      outer_failure <- function() inner_failure()
      info <- rlang::try_fetch(
        outer_failure(),
        error = capture_error_info
      )

      expect_true(any(grepl("inner_failure", info$traceback, fixed = TRUE)))
      expect_false(any(grepl("tryCatchOne", info$traceback, fixed = TRUE)))
    })

    it("preserves attached traces after the handler unwinds", {
      inner_failure <- function() stop("nested error")
      err <- rlang::try_fetch(
        inner_failure(),
        error = function(e) {
          normalized <- bayesim_fit_error(conditionMessage(e))
          normalized$trace <- rlang::trace_back()
          normalized
        }
      )

      info <- capture_error_info(err)
      expect_true(any(grepl("inner_failure", info$traceback, fixed = TRUE)))
    })
  })

  describe("flatten_with_prefix()", {
    it("flattens scalars and named vectors with a shared prefix", {
      result <- flatten_with_prefix(
        list(a = 1, b = c(x = 2, y = 3), c = 4),
        "param"
      )

      expect_equal(result$a, 1)
      expect_equal(result$param__b__x, 2)
      expect_equal(result$param__b__y, 3)
      expect_equal(result$c, 4)
    })

    it("does not flatten unnamed or single-element vectors", {
      expect_equal(
        flatten_with_prefix(list(a = c(1, 2, 3)), "test")$a,
        c(1, 2, 3)
      )
      expect_equal(
        flatten_with_prefix(list(a = c(value = 1)), "test")$a,
        c(value = 1)
      )
    })
  })

  describe("Hashing utilities", {
    it("compute_hash() is stable and discriminating", {
      expect_equal(compute_hash(list(a = 1)), compute_hash(list(a = 1)))
      expect_false(isTRUE(all.equal(
        compute_hash(list(a = 1)),
        compute_hash(list(a = 2))
      )))
    })
  })
})

# =============================================================================
# 7. Validators (from contracts.R)
# =============================================================================

describe("Validators", {
  describe("validate_data_bundle()", {
    it("accepts bundles with and without a held-out test set", {
      with_test <- list(
        train = data.frame(x = 1:10, y = rnorm(10)),
        test = data.frame(x = 1:5, y = rnorm(5)),
        response = "y",
        true_params = c(beta = 1.0, sigma = 0.5),
        vars_of_interest = c("beta", "sigma")
      )
      without_test <- list(
        train = data.frame(x = 1:10, y = rnorm(10)),
        test = NULL,
        response = "y",
        true_params = c(param1 = 1.0, param2 = 2.0),
        vars_of_interest = c("param1", "param2")
      )
      expect_silent(validate_data_bundle(with_test))
      expect_silent(validate_data_bundle(without_test))
    })

    it("rejects bundles missing train data or the response variable", {
      expect_error(validate_data_bundle(list(
        train = NULL,
        test = NULL,
        response = "y"
      )))
      expect_error(validate_data_bundle(list(
        train = data.frame(x = 1:10, y = rnorm(10)),
        test = NULL
      )))
      expect_error(validate_data_bundle("not a list"))
    })
  })

  describe("validate_fitter() / validate_metric()", {
    it("accept typed objects and reject untyped input", {
      expect_silent(validate_fitter(MockFitter()))
      expect_error(validate_fitter("not a fitter"))

      TestMetricForValidation <- S7::new_class(
        "TestMetricForValidation",
        parent = Metric,
        properties = list(
          name = S7::new_property(S7::class_character, default = "test_metric")
        )
      )
      expect_silent(validate_metric(TestMetricForValidation(name = "tm")))
      expect_error(validate_metric("not a metric"))
      expect_error(validate_metric(NULL))
    })
  })
})
