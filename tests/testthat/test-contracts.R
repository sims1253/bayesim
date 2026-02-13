# test-contracts.R
# Comprehensive tests for Phase 1 contracts and core infrastructure
# testthat 3e style with describe/it blocks

# Load the package
library(bayesim)

# =============================================================================
# 1. Error Classes (from errors.R)
# =============================================================================

describe("Error Classes", {
  describe("bayesim_error()", {
    it("creates correct class hierarchy", {
      err <- bayesim_error("Test error message")
      expect_s3_class(err, "bayesim_error")
      expect_s3_class(err, "error")
      expect_s3_class(err, "condition")
    })

    it("stores message correctly", {
      err <- bayesim_error("Test error message")
      expect_equal(err$message, "Test error message")
    })

    it("accepts optional call parameter", {
      call <- quote(my_function())
      err <- bayesim_error("Test error", call = call)
      expect_equal(err$call, call)
    })

    it("can be raised with stop()", {
      expect_error(
        stop(bayesim_error("Test error")),
        "Test error"
      )
    })
  })

  describe("is_bayesim_error()", {
    it("returns TRUE for bayesim_error objects", {
      err <- bayesim_error("Test")
      expect_true(is_bayesim_error(err))
    })

    it("returns TRUE for all derived error types", {
      expect_true(is_bayesim_error(bayesim_config_error("Test")))
      expect_true(is_bayesim_error(bayesim_contract_error("Test")))
      expect_true(is_bayesim_error(bayesim_checkpoint_error("Test")))
      expect_true(is_bayesim_error(bayesim_internal_error("Test")))
      expect_true(is_bayesim_error(bayesim_data_error("Test")))
      expect_true(is_bayesim_error(bayesim_fit_error("Test")))
      expect_true(is_bayesim_error(bayesim_metric_error("Test")))
    })

    it("returns FALSE for non-bayesim errors", {
      expect_false(is_bayesim_error(simpleError("Test")))
      expect_false(is_bayesim_error("Test"))
      expect_false(is_bayesim_error(NULL))
      expect_false(is_bayesim_error(list()))
    })
  })

  describe("Fatal errors", {
    it("bayesim_config_error() has correct class hierarchy", {
      err <- bayesim_config_error("Invalid config")
      expect_s3_class(err, "bayesim_config_error")
      expect_s3_class(err, "bayesim_error")
      expect_s3_class(err, "error")
    })

    it("bayesim_contract_error() has correct class hierarchy", {
      err <- bayesim_contract_error("Contract violation")
      expect_s3_class(err, "bayesim_contract_error")
      expect_s3_class(err, "bayesim_error")
      expect_s3_class(err, "error")
    })

    it("bayesim_checkpoint_error() has correct class hierarchy", {
      err <- bayesim_checkpoint_error("Checkpoint failed")
      expect_s3_class(err, "bayesim_checkpoint_error")
      expect_s3_class(err, "bayesim_error")
      expect_s3_class(err, "error")
    })

    it("bayesim_internal_error() has correct class hierarchy", {
      err <- bayesim_internal_error("Internal error")
      expect_s3_class(err, "bayesim_internal_error")
      expect_s3_class(err, "bayesim_error")
      expect_s3_class(err, "error")
    })
  })

  describe("Recoverable errors", {
    it("bayesim_data_error() has correct class hierarchy", {
      err <- bayesim_data_error("Data generation failed")
      expect_s3_class(err, "bayesim_data_error")
      expect_s3_class(err, "bayesim_error")
      expect_s3_class(err, "error")
    })

    it("bayesim_fit_error() has correct class hierarchy", {
      err <- bayesim_fit_error("Model fitting failed")
      expect_s3_class(err, "bayesim_fit_error")
      expect_s3_class(err, "bayesim_error")
      expect_s3_class(err, "error")
    })

    it("bayesim_metric_error() has correct class hierarchy", {
      err <- bayesim_metric_error("Metric computation failed")
      expect_s3_class(err, "bayesim_metric_error")
      expect_s3_class(err, "bayesim_error")
      expect_s3_class(err, "error")
    })
  })

  describe("is_fatal_error()", {
    it("returns TRUE for fatal error types", {
      expect_true(is_fatal_error(bayesim_config_error("Test")))
      expect_true(is_fatal_error(bayesim_contract_error("Test")))
      expect_true(is_fatal_error(bayesim_checkpoint_error("Test")))
      expect_true(is_fatal_error(bayesim_internal_error("Test")))
    })

    it("returns FALSE for recoverable error types", {
      expect_false(is_fatal_error(bayesim_data_error("Test")))
      expect_false(is_fatal_error(bayesim_fit_error("Test")))
      expect_false(is_fatal_error(bayesim_metric_error("Test")))
    })

    it("returns FALSE for base bayesim_error", {
      expect_false(is_fatal_error(bayesim_error("Test")))
    })

    it("returns FALSE for non-bayesim errors", {
      expect_false(is_fatal_error(simpleError("Test")))
      expect_false(is_fatal_error(NULL))
    })
  })

  describe("is_recoverable_error()", {
    it("returns TRUE for recoverable error types", {
      expect_true(is_recoverable_error(bayesim_data_error("Test")))
      expect_true(is_recoverable_error(bayesim_fit_error("Test")))
      expect_true(is_recoverable_error(bayesim_metric_error("Test")))
    })

    it("returns FALSE for fatal error types", {
      expect_false(is_recoverable_error(bayesim_config_error("Test")))
      expect_false(is_recoverable_error(bayesim_contract_error("Test")))
      expect_false(is_recoverable_error(bayesim_checkpoint_error("Test")))
      expect_false(is_recoverable_error(bayesim_internal_error("Test")))
    })

    it("returns FALSE for base bayesim_error", {
      expect_false(is_recoverable_error(bayesim_error("Test")))
    })

    it("returns FALSE for non-bayesim errors", {
      expect_false(is_recoverable_error(simpleError("Test")))
      expect_false(is_recoverable_error(NULL))
    })
  })

  describe("Error handling with tryCatch", {
    it("can catch bayesim_error specifically", {
      result <- tryCatch(
        stop(bayesim_config_error("Config error")),
        bayesim_error = function(e) e$message
      )
      expect_equal(result, "Config error")
    })

    it("can distinguish fatal vs recoverable in handler", {
      result <- tryCatch(
        {
          stop(bayesim_fit_error("Fit failed"))
        },
        bayesim_error = function(e) {
          if (is_fatal_error(e)) "fatal" else "recoverable"
        }
      )
      expect_equal(result, "recoverable")
    })
  })
})


# =============================================================================
# 2. Result Constructors (from results.R)
# =============================================================================

describe("Result Constructors", {
  describe("new_fit_result()", {
    it("creates valid objects with default parameters", {
      result <- new_fit_result(
        success = TRUE,
        diagnostics = list(),
        timing = list(total = 1.0)
      )
      expect_s3_class(result, "bayesim_fit_result")
      expect_true(is_bayesim_fit_result(result))
    })

    it("creates successful fit result with all fields", {
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

      expect_true(result$success)
      expect_equal(nrow(result$draws), 50)
      expect_equal(ncol(result$draws), 2)
      expect_equal(colnames(result$draws), c("alpha", "beta"))
      expect_equal(length(result$warnings), 2)
      expect_equal(result$timing$total, 10.5)
    })

    it("creates failed fit result with error", {
      err <- bayesim_fit_error("Convergence failed")
      result <- new_fit_result(
        success = FALSE,
        error = err,
        diagnostics = list(),
        timing = list(total = 2.0)
      )

      expect_false(result$success)
      expect_equal(result$error, err)
    })

    it("validates success/error consistency - error with success=TRUE fails", {
      expect_error(
        new_fit_result(
          success = TRUE,
          error = simpleError("Error"),
          diagnostics = list(),
          timing = list(total = 1.0)
        ),
        "When success is TRUE, error must be NULL"
      )
    })

    it("validates success/error consistency - no error with success=FALSE fails", {
      expect_error(
        new_fit_result(
          success = FALSE,
          error = NULL,
          diagnostics = list(),
          timing = list(total = 1.0)
        ),
        "When success is FALSE, error must be non-NULL"
      )
    })

    it("validates timing is a list", {
      # Note: new_fit_result tries to access timing$total before validation,
      # so we get "$ operator is invalid for atomic vectors" error
      expect_error(
        new_fit_result(
          success = TRUE,
          diagnostics = list(),
          timing = "not a list"
        ),
        "\\$ operator is invalid"
      )
    })

    it("validates timing$total is scalar numeric", {
      expect_error(
        new_fit_result(
          success = TRUE,
          diagnostics = list(),
          timing = list(total = c(1, 2))
        ),
        "timing\\$total must be a scalar numeric"
      )
    })

    it("validates timing$total is non-negative", {
      expect_error(
        new_fit_result(
          success = TRUE,
          diagnostics = list(),
          timing = list(total = -1.0)
        ),
        "timing\\$total must be >= 0"
      )
    })

    it("validates draws is a matrix with colnames", {
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      # No colnames - should fail
      expect_error(
        new_fit_result(
          success = TRUE,
          draws = draws,
          diagnostics = list(),
          timing = list(total = 1.0)
        ),
        "draws matrix must have colnames"
      )
    })

    it("validates draws must be a matrix", {
      expect_error(
        new_fit_result(
          success = TRUE,
          draws = data.frame(a = 1:10),
          diagnostics = list(),
          timing = list(total = 1.0)
        ),
        "draws must be a matrix when not NULL"
      )
    })

    it("validates warnings is character vector", {
      expect_error(
        new_fit_result(
          success = TRUE,
          diagnostics = list(),
          timing = list(total = 1.0),
          warnings = list("not", "character")
        ),
        "warnings must be a character vector"
      )
    })

    it("validates diagnostics is a list", {
      expect_error(
        new_fit_result(
          success = TRUE,
          diagnostics = "not a list",
          timing = list(total = 1.0)
        ),
        "diagnostics must be a list"
      )
    })

    it("handles NULL warnings by converting to character()", {
      result <- new_fit_result(
        success = TRUE,
        diagnostics = list(),
        timing = list(total = 1.0),
        warnings = NULL
      )
      expect_equal(result$warnings, character())
    })

    it("handles NULL timing by creating default", {
      result <- new_fit_result(
        success = TRUE,
        diagnostics = list(),
        timing = NULL
      )
      expect_equal(result$timing$total, 0)
    })
  })

  describe("new_task_result()", {
    it("creates valid success task result", {
      result <- new_task_result(
        task_id = "d001_f001_r00001",
        status = "success",
        metrics = list(rmse = 0.05, bias = 0.01),
        diagnostics = list(n_eff = 500),
        timing = list(total = 5.2)
      )

      expect_s3_class(result, "bayesim_task_result")
      expect_true(is_bayesim_task_result(result))
      expect_equal(result$task_id, "d001_f001_r00001")
      expect_equal(result$status, "success")
    })

    it("creates valid failed task result", {
      result <- new_task_result(
        task_id = "d001_f001_r00002",
        status = "failed",
        error = list(
          error_class = "convergence_error",
          error_message = "R-hat > 1.1"
        ),
        timing = list(total = 2.0)
      )

      expect_equal(result$status, "failed")
      expect_equal(result$error$error_class, "convergence_error")
    })

    it("creates valid skipped task result", {
      result <- new_task_result(
        task_id = "d001_f001_r00003",
        status = "skipped",
        timing = list(total = 0.0)
      )

      expect_equal(result$status, "skipped")
    })

    it("validates task_id is scalar character", {
      expect_error(
        new_task_result(
          task_id = 123,
          status = "success",
          metrics = list(a = 1),
          timing = list(total = 1.0)
        ),
        "task_id must be a scalar character"
      )
    })

    it("validates status is valid value", {
      expect_error(
        new_task_result(
          task_id = "test",
          status = "invalid",
          timing = list(total = 1.0)
        ),
        "status must be one of"
      )
    })

    it("validates success requires non-NULL metrics", {
      expect_error(
        new_task_result(
          task_id = "test",
          status = "success",
          metrics = NULL,
          timing = list(total = 1.0)
        ),
        "When status is 'success', metrics must not be NULL"
      )
    })

    it("validates failed requires non-NULL error", {
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

    it("validates timing$total is non-negative", {
      expect_error(
        new_task_result(
          task_id = "test",
          status = "skipped",
          timing = list(total = -1.0)
        ),
        "timing\\$total must be >= 0"
      )
    })

    it("validates warnings is character vector", {
      expect_error(
        new_task_result(
          task_id = "test",
          status = "skipped",
          timing = list(total = 1.0),
          warnings = list(a = 1)
        ),
        "warnings must be a character vector"
      )
    })
  })

  describe("new_simulation_result()", {
    it("creates valid simulation result", {
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
      expect_true(is_bayesim_simulation_result(result))
      expect_equal(result$config_fingerprint, "abc123")
      expect_equal(length(result$task_results), 1)
    })

    it("validates config_fingerprint is scalar character", {
      expect_error(
        new_simulation_result(
          config_fingerprint = c("a", "b"),
          task_results = list()
        ),
        "config_fingerprint must be a scalar character"
      )
    })

    it("validates task_results is a list", {
      expect_error(
        new_simulation_result(
          config_fingerprint = "test",
          task_results = "not a list"
        ),
        "task_results must be a list"
      )
    })

    it("validates all task_results elements are bayesim_task_result", {
      expect_error(
        new_simulation_result(
          config_fingerprint = "test",
          task_results = list(list(not_a_result = TRUE))
        ),
        "All elements of task_results must be bayesim_task_result objects"
      )
    })

    it("validates summary is a data.frame", {
      task <- new_task_result(
        task_id = "test",
        status = "skipped",
        timing = list(total = 0)
      )
      expect_error(
        new_simulation_result(
          config_fingerprint = "test",
          task_results = list(task),
          summary = "not a data.frame"
        ),
        "summary must be a data.frame or tibble"
      )
    })

    it("validates errors is a data.frame", {
      task <- new_task_result(
        task_id = "test",
        status = "skipped",
        timing = list(total = 0)
      )
      expect_error(
        new_simulation_result(
          config_fingerprint = "test",
          task_results = list(task),
          summary = data.frame(),
          errors = "not a data.frame"
        ),
        "errors must be a data.frame or tibble"
      )
    })

    it("validates checkpoint_path is NULL or scalar character", {
      task <- new_task_result(
        task_id = "test",
        status = "skipped",
        timing = list(total = 0)
      )
      expect_error(
        new_simulation_result(
          config_fingerprint = "test",
          task_results = list(task),
          checkpoint_path = c("a", "b")
        ),
        "checkpoint_path must be NULL or a scalar character"
      )
    })

    it("accepts NULL checkpoint_path", {
      task <- new_task_result(
        task_id = "test",
        status = "skipped",
        timing = list(total = 0)
      )
      result <- new_simulation_result(
        config_fingerprint = "test",
        task_results = list(task),
        checkpoint_path = NULL
      )
      expect_null(result$checkpoint_path)
    })

    it("handles NULL inputs with defaults", {
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

  describe("Print methods", {
    it("print.bayesim_fit_result() works without error", {
      result <- new_fit_result(
        success = TRUE,
        diagnostics = list(rhat = 1.01),
        timing = list(total = 5.5)
      )
      expect_output(print(result), "<bayesim_fit_result>")
    })

    it("print.bayesim_task_result() works without error", {
      result <- new_task_result(
        task_id = "test",
        status = "success",
        metrics = list(rmse = 0.05),
        timing = list(total = 5.0)
      )
      expect_output(print(result), "<bayesim_task_result>")
    })

    it("print.bayesim_simulation_result() works without error", {
      task <- new_task_result(
        task_id = "test",
        status = "success",
        metrics = list(rmse = 0.05),
        timing = list(total = 5.0)
      )
      result <- new_simulation_result(
        config_fingerprint = "test",
        task_results = list(task)
      )
      expect_output(print(result), "<bayesim_simulation_result>")
    })
  })
})


# =============================================================================
# 3. Fitter Class (from fitter.R)
# =============================================================================

describe("Fitter Class", {
  describe("Fitter abstract class", {
    it("is abstract", {
      expect_true(Fitter@abstract)
    })

    it("has expected properties", {
      expect_true("name" %in% names(Fitter@properties))
      expect_true("supports_predictions" %in% names(Fitter@properties))
      expect_true("supports_log_lik" %in% names(Fitter@properties))
      expect_true("supports_loo" %in% names(Fitter@properties))
    })

    it("cannot be instantiated directly", {
      expect_error(
        Fitter(name = "test"),
        "abstract"
      )
    })
  })

  describe("MockFitter", {
    # Create fitter in each test instead of using before_each
    it("inherits from Fitter", {
      fitter <- MockFitter()
      expect_true(S7::S7_inherits(fitter))
      expect_true(S7::S7_inherits(fitter, Fitter))
    })

    it("has default name 'mock'", {
      fitter <- MockFitter()
      expect_equal(fitter@name, "mock")
    })

    it("fit() returns bayesim_fit_result", {
      fitter <- MockFitter()
      data_bundle <- list(
        train = data.frame(x = 1:10, y = rnorm(10)),
        test = NULL,
        response = "y"
      )
      fit_spec <- data.frame(model = "test")
      task_ctx <- list(task_id = "test_task")

      result <- fit(
        fitter,
        data_bundle,
        fit_spec,
        seed = 123L,
        task_ctx = task_ctx
      )

      expect_s3_class(result, "bayesim_fit_result")
      expect_equal(result$fit$fitter, "mock")
    })

    it("extract_draws() returns matrix with colnames", {
      fitter <- MockFitter()
      # First create a fit result
      data_bundle <- list(
        train = data.frame(x = 1:10, y = rnorm(10)),
        test = NULL,
        response = "y"
      )
      fit_result <- fit(
        fitter,
        data_bundle,
        fit_spec = data.frame(model = "test"),
        seed = 42L,
        task_ctx = list(task_id = "test")
      )

      draws <- extract_draws(fitter, fit_result)

      expect_true(is.matrix(draws))
      expect_false(is.null(colnames(draws)))
      expect_true(all(c("intercept", "beta", "sigma") %in% colnames(draws)))
    })

    it("extract_draws() respects variables argument", {
      fitter <- MockFitter()
      data_bundle <- list(
        train = data.frame(x = 1:10, y = rnorm(10)),
        test = NULL,
        response = "y"
      )
      fit_result <- fit(
        fitter,
        data_bundle,
        fit_spec = data.frame(model = "test"),
        seed = 42L,
        task_ctx = list(task_id = "test")
      )

      draws <- extract_draws(fitter, fit_result, variables = c("beta", "sigma"))

      expect_equal(colnames(draws), c("beta", "sigma"))
    })

    it("predict_fit() returns proper structure", {
      fitter <- MockFitter()
      data_bundle <- list(
        train = data.frame(x = 1:10, y = rnorm(10)),
        test = NULL,
        response = "y"
      )
      fit_result <- fit(
        fitter,
        data_bundle,
        fit_spec = data.frame(model = "test"),
        seed = 42L,
        task_ctx = list(task_id = "test")
      )

      preds <- predict_fit(fitter, fit_result)

      expect_true(is.list(preds))
      expect_true("predicted_mean" %in% names(preds))
      expect_true("predicted_samples" %in% names(preds))
      expect_true("predicted_sd" %in% names(preds))
      expect_equal(length(preds$predicted_mean), 10)
    })

    it("log_lik() returns matrix", {
      fitter <- MockFitter()
      data_bundle <- list(
        train = data.frame(x = 1:10, y = rnorm(10)),
        test = NULL,
        response = "y"
      )
      fit_result <- fit(
        fitter,
        data_bundle,
        fit_spec = data.frame(model = "test"),
        seed = 42L,
        task_ctx = list(task_id = "test")
      )

      ll <- log_lik(fitter, fit_result)

      expect_true(is.matrix(ll))
      expect_equal(nrow(ll), 10) # 10 observations
    })

    it("loo() returns proper structure", {
      fitter <- MockFitter()
      data_bundle <- list(
        train = data.frame(x = 1:10, y = rnorm(10)),
        test = NULL,
        response = "y"
      )
      fit_result <- fit(
        fitter,
        data_bundle,
        fit_spec = data.frame(model = "test"),
        seed = 42L,
        task_ctx = list(task_id = "test")
      )

      loo_result <- loo(fitter, fit_result)

      expect_true(is.list(loo_result))
      expect_true("elpd" %in% names(loo_result))
      expect_true("p_loo" %in% names(loo_result))
      expect_true("elpd_se" %in% names(loo_result))
      expect_true("pareto_k" %in% names(loo_result))
    })

    it("diagnostics() returns named list", {
      fitter <- MockFitter()
      data_bundle <- list(
        train = data.frame(x = 1:10, y = rnorm(10)),
        test = NULL,
        response = "y"
      )
      fit_result <- fit(
        fitter,
        data_bundle,
        fit_spec = data.frame(model = "test"),
        seed = 42L,
        task_ctx = list(task_id = "test")
      )

      diag <- diagnostics(fitter, fit_result)

      expect_true(is.list(diag))
      expect_true("rhat_max" %in% names(diag))
      expect_true("ess_bulk" %in% names(diag))
      expect_true("divergent" %in% names(diag))
    })
  })

  describe("Custom Fitter implementation", {
    # Define TestFitter inline for each test
    it("can extend Fitter with custom implementation", {
      TestFitter <- S7::new_class(
        "TestFitter",
        parent = Fitter,
        properties = list(
          name = S7::new_property(S7::class_character, default = "test")
        )
      )
      # Register method separately using S7 generics-based system
      S7::method(fit, TestFitter) <- function(
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
          diagnostics = list(),
          timing = list(total = 1.0)
        )
      }

      fitter <- TestFitter()
      expect_true(S7::S7_inherits(fitter))
      expect_equal(fitter@name, "test")
    })

    it("unimplemented methods throw S7 method not found error", {
      TestFitter <- S7::new_class(
        "TestFitter",
        parent = Fitter,
        properties = list(
          name = S7::new_property(S7::class_character, default = "test")
        )
      )
      S7::method(fit, TestFitter) <- function(
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
          diagnostics = list(),
          timing = list(total = 1.0)
        )
      }

      fitter <- TestFitter()
      fit_result <- fit(
        fitter,
        data.frame(x = 1),
        data.frame(),
        1L,
        list(task_id = "test")
      )

      # S7 throws "Can't find method for" error when method is not implemented
      expect_error(
        extract_draws(fitter, fit_result),
        "Can't find method"
      )
    })
  })
})


# =============================================================================
# 4. Metric Class (from metric.R)
# =============================================================================

describe("Metric Class", {
  describe("Metric abstract class", {
    it("has expected properties", {
      expect_true("name" %in% names(Metric@properties))
      expect_true("needs" %in% names(Metric@properties))
      expect_true("required" %in% names(Metric@properties))
    })

    it("compute() throws S7 method not found error for unimplemented methods", {
      TestMetric <- S7::new_class(
        "TestMetric",
        parent = Metric,
        properties = list(
          name = S7::new_property(S7::class_character, default = "test")
        )
      )
      metric <- TestMetric()

      # S7 throws "Can't find method for" error when method is not implemented
      expect_error(
        compute(metric, NULL, NULL, NULL, NULL),
        "Can't find method"
      )
    })
  })

  describe("validate_metric_output()", {
    it("accepts valid scalar outputs", {
      output <- list(rmse = 0.5, n_obs = 100L, success = TRUE)
      expect_silent(validate_metric_output(output, "test_metric"))
    })

    it("accepts named numeric vectors", {
      output <- list(params = c(alpha = 0.1, beta = 0.2, gamma = 0.3))
      expect_silent(validate_metric_output(output, "test_metric"))
    })

    it("accepts mixed scalar and named vector outputs", {
      output <- list(
        rmse = 0.5,
        params = c(alpha = 0.1, beta = 0.2)
      )
      expect_silent(validate_metric_output(output, "test_metric"))
    })

    it("rejects non-list output", {
      expect_error(
        validate_metric_output(c(a = 1, b = 2), "test_metric"),
        "must be a list"
      )
    })

    it("rejects empty list", {
      expect_error(
        validate_metric_output(list(), "test_metric"),
        "cannot be empty"
      )
    })

    it("rejects unnamed elements", {
      expect_error(
        validate_metric_output(list(0.5, 1.0), "test_metric"),
        "unnamed or empty-named"
      )
    })

    it("rejects elements with empty names", {
      output <- list(rmse = 0.5)
      names(output) <- ""
      expect_error(
        validate_metric_output(output, "test_metric"),
        "unnamed or empty-named"
      )
    })

    it("rejects NULL values", {
      expect_error(
        validate_metric_output(list(rmse = NULL), "test_metric"),
        "is NULL \\(not allowed\\)"
      )
    })

    it("rejects nested lists", {
      expect_error(
        validate_metric_output(list(nested = list(a = 1)), "test_metric"),
        "must be either"
      )
    })

    it("rejects data frames", {
      expect_error(
        validate_metric_output(list(df = data.frame(a = 1:3)), "test_metric"),
        "must be either"
      )
    })

    it("rejects matrices", {
      expect_error(
        validate_metric_output(list(mat = matrix(1:4, 2, 2)), "test_metric"),
        "must be either"
      )
    })

    it("rejects unnamed numeric vectors", {
      expect_error(
        validate_metric_output(list(params = c(1, 2, 3)), "test_metric"),
        "unnamed numeric vector"
      )
    })

    it("rejects vectors of other types", {
      expect_error(
        validate_metric_output(c(TRUE, FALSE, TRUE), "test_metric"),
        "must be a list"
      )
    })
  })

  describe("flatten_metric_output()", {
    it("produces correct names for scalar values", {
      output <- list(rmse = 0.5, n_obs = 100L)
      flat <- flatten_metric_output(output, "my_metric")

      expect_equal(flat$my_metric__rmse, 0.5)
      expect_equal(flat$my_metric__n_obs, 100L)
    })

    it("handles named vectors correctly", {
      output <- list(params = c(alpha = 0.1, beta = 0.2, gamma = 0.3))
      flat <- flatten_metric_output(output, "estimates")

      expect_equal(flat$estimates__params__alpha, 0.1)
      expect_equal(flat$estimates__params__beta, 0.2)
      expect_equal(flat$estimates__params__gamma, 0.3)
    })

    it("handles mixed scalar and named vectors", {
      output <- list(
        rmse = 0.5,
        params = c(alpha = 0.1, beta = 0.2)
      )
      flat <- flatten_metric_output(output, "test")

      expect_equal(flat$test__rmse, 0.5)
      expect_equal(flat$test__params__alpha, 0.1)
      expect_equal(flat$test__params__beta, 0.2)
    })

    it("validates input before flattening", {
      expect_error(
        flatten_metric_output(list(unnamed = NULL), "test"),
        "is NULL"
      )
    })
  })

  describe("Custom Metric implementation", {
    # Define RMSEMetric inline for each test
    it("can extend Metric with custom implementation", {
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
      # Register method separately using S7 generics-based system
      S7::method(compute, RMSEMetric) <- function(
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
      expect_true(S7::S7_inherits(metric))
      # S7 default values for class_character with default = "string" work correctly
      expect_true(is.character(metric@name))
    })

    it("compute() returns valid output", {
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
      S7::method(compute, RMSEMetric) <- function(
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
      fit_result <- list()
      data_bundle <- list(
        test = data.frame(y = c(1, 2, 3)),
        response = "y"
      )
      context <- list(predictions = list(predicted_mean = c(1.1, 1.9, 3.2)))
      task_ctx <- list(task_id = "test")

      output <- compute(metric, fit_result, data_bundle, context, task_ctx)

      expect_silent(validate_metric_output(output, "rmse"))
      expect_true("value" %in% names(output))
      expect_true("n_obs" %in% names(output))
    })
  })
})


# =============================================================================
# 5. SimulationConfig (from simulation-config.R)
# =============================================================================

describe("SimulationConfig", {
  test_data_gen <- function(data_spec, seed, task_ctx) {
    list(
      train = data.frame(x = 1:10, y = rnorm(10)),
      test = NULL,
      response = "y"
    )
  }

  describe("simulation_config()", {
    it("creates valid objects with required parameters", {
      config <- simulation_config(
        data_grid = data.frame(n = 100),
        fit_grid = data.frame(model = "baseline"),
        data_generator = test_data_gen,
        seed = 42L
      )

      expect_true(S7::S7_inherits(config))
      expect_true(is_simulation_config(config))
    })

    it("creates valid objects with all parameters", {
      config <- simulation_config(
        data_grid = data.frame(n = c(100, 500)),
        fit_grid = data.frame(model = c("A", "B")),
        data_generator = test_data_gen,
        fitter = MockFitter(),
        metrics = list(),
        n_replicates = 100L,
        seed = 42L,
        result_path = "/tmp/results",
        checkpoint_every = 25L,
        retain = c("metrics", "diagnostics"),
        max_errors = 10
      )

      expect_true(S7::S7_inherits(config))
      expect_equal(config@n_replicates, 100L)
      expect_equal(config@seed, 42L)
      expect_equal(config@result_path, "/tmp/results")
      expect_equal(config@checkpoint_every, 25L)
      expect_equal(config@retain, c("metrics", "diagnostics"))
    })

    it("validates data_grid is a data.frame", {
      expect_error(
        simulation_config(
          data_grid = "not a data.frame",
          fit_grid = data.frame(),
          data_generator = test_data_gen,
          seed = 42L
        ),
        "data_grid must be a data.frame"
      )
    })

    it("validates data_grid has at least 1 row", {
      expect_error(
        simulation_config(
          data_grid = data.frame(),
          fit_grid = data.frame(a = 1),
          data_generator = test_data_gen,
          seed = 42L
        ),
        "data_grid must have at least 1 row"
      )
    })

    it("validates fit_grid is a data.frame", {
      expect_error(
        simulation_config(
          data_grid = data.frame(a = 1),
          fit_grid = "not a data.frame",
          data_generator = test_data_gen,
          seed = 42L
        ),
        "fit_grid must be a data.frame"
      )
    })

    it("validates fit_grid has at least 1 row", {
      expect_error(
        simulation_config(
          data_grid = data.frame(a = 1),
          fit_grid = data.frame(),
          data_generator = test_data_gen,
          seed = 42L
        ),
        "fit_grid must have at least 1 row"
      )
    })

    it("validates data_generator is a function", {
      expect_error(
        simulation_config(
          data_grid = data.frame(a = 1),
          fit_grid = data.frame(a = 1),
          data_generator = "not a function",
          seed = 42L
        ),
        "data_generator must be a function"
      )
    })

    it("validates data_generator has at least 3 arguments", {
      bad_gen <- function(a, b) list()
      expect_error(
        simulation_config(
          data_grid = data.frame(a = 1),
          fit_grid = data.frame(a = 1),
          data_generator = bad_gen,
          seed = 42L
        ),
        "data_generator must accept at least 3 arguments"
      )
    })

    it("validates fitter is S7 object", {
      expect_error(
        simulation_config(
          data_grid = data.frame(a = 1),
          fit_grid = data.frame(a = 1),
          data_generator = test_data_gen,
          fitter = "not a fitter",
          seed = 42L
        ),
        "fitter must be an S7 Fitter object"
      )
    })

    it("validates n_replicates is positive integer", {
      expect_error(
        simulation_config(
          data_grid = data.frame(a = 1),
          fit_grid = data.frame(a = 1),
          data_generator = test_data_gen,
          n_replicates = 0L,
          seed = 42L
        ),
        "n_replicates must be a positive integer"
      )
    })

    it("validates seed is single integer", {
      expect_error(
        simulation_config(
          data_grid = data.frame(a = 1),
          fit_grid = data.frame(a = 1),
          data_generator = test_data_gen,
          seed = c(1, 2)
        ),
        "seed must be a single integer"
      )
    })

    it("validates result_path is NULL or single character", {
      expect_error(
        simulation_config(
          data_grid = data.frame(a = 1),
          fit_grid = data.frame(a = 1),
          data_generator = test_data_gen,
          result_path = c("a", "b"),
          seed = 42L
        ),
        "result_path must be NULL or a single character string"
      )
    })

    it("validates checkpoint_every is positive integer", {
      expect_error(
        simulation_config(
          data_grid = data.frame(a = 1),
          fit_grid = data.frame(a = 1),
          data_generator = test_data_gen,
          checkpoint_every = 0L,
          seed = 42L
        ),
        "checkpoint_every must be a positive integer"
      )
    })

    it("validates retain contains valid options", {
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
    })

    it("validates max_errors is single numeric", {
      expect_error(
        simulation_config(
          data_grid = data.frame(a = 1),
          fit_grid = data.frame(a = 1),
          data_generator = test_data_gen,
          max_errors = c(1, 2),
          seed = 42L
        ),
        "max_errors must be a single numeric value"
      )
    })

    it("validates max_errors is non-negative or Inf", {
      expect_error(
        simulation_config(
          data_grid = data.frame(a = 1),
          fit_grid = data.frame(a = 1),
          data_generator = test_data_gen,
          max_errors = -5,
          seed = 42L
        ),
        "max_errors must be Inf or a non-negative number"
      )
    })

    it("accepts Inf for max_errors", {
      config <- simulation_config(
        data_grid = data.frame(a = 1),
        fit_grid = data.frame(a = 1),
        data_generator = test_data_gen,
        max_errors = Inf,
        seed = 42L
      )
      expect_equal(config@max_errors, Inf)
    })
  })

  describe("is_simulation_config()", {
    it("returns TRUE for SimulationConfig objects", {
      config <- simulation_config(
        data_grid = data.frame(a = 1),
        fit_grid = data.frame(b = 1),
        data_generator = test_data_gen,
        seed = 42L
      )
      expect_true(is_simulation_config(config))
    })

    it("returns FALSE for non-SimulationConfig objects", {
      expect_false(is_simulation_config(list()))
      expect_false(is_simulation_config(NULL))
      expect_false(is_simulation_config("string"))
    })
  })

  describe("as_config_spec()", {
    it("produces correct structure", {
      config <- simulation_config(
        data_grid = data.frame(n = 100),
        fit_grid = data.frame(model = "baseline"),
        data_generator = test_data_gen,
        n_replicates = 10L,
        seed = 42L
      )

      spec <- as_config_spec(config)

      expect_true(is.list(spec))
      expect_true("data_grid" %in% names(spec))
      expect_true("fit_grid" %in% names(spec))
      expect_true("n_replicates" %in% names(spec))
      expect_true("seed" %in% names(spec))
    })

    it("excludes runtime-specific fields", {
      config <- simulation_config(
        data_grid = data.frame(a = 1),
        fit_grid = data.frame(b = 1),
        data_generator = test_data_gen,
        result_path = "/tmp/test",
        checkpoint_every = 100L,
        seed = 42L
      )

      spec <- as_config_spec(config)

      expect_false("result_path" %in% names(spec))
      expect_false("checkpoint_every" %in% names(spec))
    })

    it("errors on non-SimulationConfig input", {
      expect_error(
        as_config_spec(list()),
        "config must be a SimulationConfig object"
      )
    })
  })

  describe("compute_config_fingerprint()", {
    it("produces consistent hashes for same config", {
      config <- simulation_config(
        data_grid = data.frame(n = 100),
        fit_grid = data.frame(model = "baseline"),
        data_generator = test_data_gen,
        n_replicates = 10L,
        seed = 42L
      )

      fp1 <- compute_config_fingerprint(config)
      fp2 <- compute_config_fingerprint(config)

      expect_equal(fp1, fp2)
      expect_true(is.character(fp1))
      expect_equal(nchar(fp1), 64) # SHA256 produces 64 hex characters
    })

    it("produces different hashes for different configs", {
      config1 <- simulation_config(
        data_grid = data.frame(n = 100),
        fit_grid = data.frame(model = "A"),
        data_generator = test_data_gen,
        seed = 42L
      )
      config2 <- simulation_config(
        data_grid = data.frame(n = 200),
        fit_grid = data.frame(model = "A"),
        data_generator = test_data_gen,
        seed = 42L
      )

      fp1 <- compute_config_fingerprint(config1)
      fp2 <- compute_config_fingerprint(config2)

      expect_true(!identical(fp1, fp2) && !isTRUE(all.equal(fp1, fp2)))
    })

    it("errors on non-SimulationConfig input", {
      expect_error(
        compute_config_fingerprint(list()),
        "config must be a SimulationConfig object"
      )
    })
  })

  describe("get_total_tasks()", {
    it("calculates correctly for simple config", {
      config <- simulation_config(
        data_grid = data.frame(n = c(100, 500)), # 2 rows
        fit_grid = data.frame(model = c("A", "B")), # 2 rows
        data_generator = test_data_gen,
        n_replicates = 100L,
        seed = 42L
      )

      expect_equal(get_total_tasks(config), 2 * 2 * 100) # 400
    })

    it("calculates correctly for single row grids", {
      config <- simulation_config(
        data_grid = data.frame(n = 100), # 1 row
        fit_grid = data.frame(model = "A"), # 1 row
        data_generator = test_data_gen,
        n_replicates = 50L,
        seed = 42L
      )

      expect_equal(get_total_tasks(config), 1 * 1 * 50) # 50
    })

    it("errors on non-SimulationConfig input", {
      expect_error(
        get_total_tasks(list()),
        "config must be a SimulationConfig object"
      )
    })
  })
})


# =============================================================================
# 6. Utility Functions (from utils.R)
# =============================================================================

describe("Utility Functions", {
  describe("format_task_id() and parse_task_id()", {
    it("format_task_id() creates correct format", {
      expect_equal(format_task_id(1, 2, 100), "d001_f002_r00100")
      expect_equal(format_task_id(999, 999, 99999), "d999_f999_r99999")
    })

    it("parse_task_id() extracts indices correctly", {
      result <- parse_task_id("d001_f002_r00100")
      expect_equal(result$data_idx, 1L)
      expect_equal(result$fit_idx, 2L)
      expect_equal(result$rep_idx, 100L)
    })

    it("format_task_id() and parse_task_id() are inverses", {
      test_cases <- list(
        c(1, 1, 1),
        c(5, 10, 100),
        c(123, 456, 7890),
        c(999, 999, 99999)
      )

      for (case in test_cases) {
        task_id <- format_task_id(case[1], case[2], case[3])
        parsed <- parse_task_id(task_id)

        expect_equal(parsed$data_idx, as.integer(case[1]))
        expect_equal(parsed$fit_idx, as.integer(case[2]))
        expect_equal(parsed$rep_idx, as.integer(case[3]))
      }
    })
  })

  describe("make_timer()", {
    it("tracks time correctly", {
      timer <- make_timer()

      # Before starting
      expect_equal(timer$elapsed(), 0)

      timer$start()
      Sys.sleep(0.1)
      timer$stop()

      elapsed <- timer$elapsed()
      expect_gte(elapsed, 0.1)
      expect_lt(elapsed, 0.5) # Should not take too long
    })

    it("returns 0 when not started", {
      timer <- make_timer()
      expect_equal(timer$elapsed(), 0)
    })

    it("returns current time when not stopped", {
      timer <- make_timer()
      timer$start()
      Sys.sleep(0.05)

      elapsed <- timer$elapsed()
      expect_gte(elapsed, 0.05)
    })

    it("can be restarted", {
      timer <- make_timer()
      timer$start()
      Sys.sleep(0.05)
      timer$stop()

      first_elapsed <- timer$elapsed()

      timer$start()
      Sys.sleep(0.05)
      timer$stop()

      second_elapsed <- timer$elapsed()
      expect_lt(second_elapsed, first_elapsed + 0.1)
    })
  })

  describe("capture_error_info()", {
    it("captures error details", {
      err <- tryCatch(
        stop("Test error message"),
        error = function(e) e
      )

      info <- capture_error_info(err)

      expect_true(is.list(info))
      expect_true("class" %in% names(info))
      expect_true("message" %in% names(info))
      expect_true("call" %in% names(info))
      expect_true("traceback" %in% names(info))
      expect_equal(info$message, "Test error message")
    })

    it("captures error class", {
      err <- tryCatch(
        stop(bayesim_config_error("Config error")),
        error = function(e) e
      )

      info <- capture_error_info(err)
      expect_true(grepl("bayesim_config_error", info$class))
    })

    it("handles errors without call", {
      err <- simpleError("Error without call")
      info <- capture_error_info(err)

      expect_equal(info$message, "Error without call")
      expect_null(info$call)
    })

    it("limits traceback length", {
      # Create a deeply nested call to generate long traceback
      err <- tryCatch(
        {
          f1 <- function() f2()
          f2 <- function() f3()
          f3 <- function() stop("Nested error")
          f1()
        },
        error = function(e) e
      )

      info <- capture_error_info(err)
      # traceback should be limited
      expect_lte(length(info$traceback), 20)
    })
  })

  describe("flatten_with_prefix()", {
    it("flattens simple list without nested vectors", {
      x <- list(a = 1, b = 2, c = 3)
      result <- flatten_with_prefix(x, "test")

      expect_equal(result$a, 1)
      expect_equal(result$b, 2)
      expect_equal(result$c, 3)
    })

    it("flattens named numeric vectors with prefix", {
      x <- list(b = c(x = 2, y = 3))
      result <- flatten_with_prefix(x, "param")

      # flatten_with_prefix uses "__" separator for all parts
      expect_equal(result$param__b__x, 2)
      expect_equal(result$param__b__y, 3)
    })

    it("handles mixed list correctly", {
      x <- list(a = 1, b = c(x = 2, y = 3), c = 4)
      result <- flatten_with_prefix(x, "param")

      expect_equal(result$a, 1)
      # flatten_with_prefix uses "__" separator for all parts
      expect_equal(result$param__b__x, 2)
      expect_equal(result$param__b__y, 3)
      expect_equal(result$c, 4)
    })

    it("does not flatten unnamed vectors", {
      x <- list(a = c(1, 2, 3)) # unnamed vector
      result <- flatten_with_prefix(x, "test")

      expect_equal(result$a, c(1, 2, 3))
    })

    it("does not flatten single-element vectors", {
      x <- list(a = c(value = 1))
      result <- flatten_with_prefix(x, "test")

      expect_equal(result$a, c(value = 1))
    })
  })

  describe("Hashing utilities", {
    it("compute_hash() produces consistent hashes", {
      obj <- list(a = 1, b = "test")
      h1 <- compute_hash(obj)
      h2 <- compute_hash(obj)

      expect_equal(h1, h2)
      expect_true(is.character(h1))
    })

    it("compute_hash() produces different hashes for different objects", {
      h1 <- compute_hash(list(a = 1))
      h2 <- compute_hash(list(a = 2))

      expect_true(!identical(h1, h2) && !isTRUE(all.equal(h1, h2)))
    })
  })
})


# =============================================================================
# 7. Validators (from contracts.R)
# =============================================================================

describe("Validators", {
  describe("validate_data_bundle()", {
    # Check if function exists
    validator_exists <- exists("validate_data_bundle", mode = "function")

    it("accepts valid bundles", {
      skip_if_not(validator_exists, "validate_data_bundle not implemented")

      bundle <- list(
        train = data.frame(x = 1:10, y = rnorm(10)),
        test = data.frame(x = 1:5, y = rnorm(5)),
        response = "y",
        true_params = c(beta = 1.0, sigma = 0.5),
        vars_of_interest = c("beta", "sigma"),
        references = c(beta = 0, sigma = 1)
      )
      expect_silent(validate_data_bundle(bundle))
    })

    it("accepts bundles without test data", {
      skip_if_not(validator_exists, "validate_data_bundle not implemented")

      bundle <- list(
        train = data.frame(x = 1:10, y = rnorm(10)),
        test = NULL,
        response = "y",
        true_params = c(param1 = 1.0, param2 = 2.0),
        vars_of_interest = c("param1", "param2"),
        references = c(param1 = 0.0, param2 = 0.0)
      )
      expect_silent(validate_data_bundle(bundle))
    })

    it("rejects bundles without train data", {
      skip_if_not(validator_exists, "validate_data_bundle not implemented")

      bundle <- list(
        train = NULL,
        test = NULL,
        response = "y"
      )
      expect_error(validate_data_bundle(bundle))
    })

    it("rejects bundles without response", {
      skip_if_not(validator_exists, "validate_data_bundle not implemented")

      bundle <- list(
        train = data.frame(x = 1:10, y = rnorm(10)),
        test = NULL
      )
      expect_error(validate_data_bundle(bundle))
    })

    it("rejects non-list input", {
      skip_if_not(validator_exists, "validate_data_bundle not implemented")
      expect_error(validate_data_bundle("not a list"))
    })
  })

  describe("validate_fitter_interface()", {
    validator_exists <- exists("validate_fitter_interface", mode = "function")

    it("accepts valid Fitter objects", {
      skip_if_not(validator_exists, "validate_fitter_interface not implemented")

      fitter <- MockFitter()
      expect_silent(validate_fitter_interface(fitter))
    })

    it("rejects non-Fitter objects", {
      skip_if_not(validator_exists, "validate_fitter_interface not implemented")
      expect_error(validate_fitter_interface("not a fitter"))
    })

    it("rejects NULL", {
      skip_if_not(validator_exists, "validate_fitter_interface not implemented")
      expect_error(validate_fitter_interface(NULL))
    })
  })

  describe("validate_metric_interface()", {
    validator_exists <- exists("validate_metric_interface", mode = "function")

    it("accepts valid Metric objects", {
      skip_if_not(validator_exists, "validate_metric_interface not implemented")

      # Define a test metric inline with name property set
      TestMetricForValidation <- S7::new_class(
        "TestMetricForValidation",
        parent = Metric,
        properties = list(
          name = S7::new_property(S7::class_character, default = "test_metric")
        )
      )
      metric <- TestMetricForValidation(name = "test_metric")
      expect_silent(validate_metric_interface(metric))
    })

    it("rejects non-Metric objects", {
      skip_if_not(validator_exists, "validate_metric_interface not implemented")
      expect_error(validate_metric_interface("not a metric"))
    })

    it("rejects NULL", {
      skip_if_not(validator_exists, "validate_metric_interface not implemented")
      expect_error(validate_metric_interface(NULL))
    })
  })
})
