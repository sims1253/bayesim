# test-regression.R
# Regression tests for issues found during code review
# testthat 3e style with describe/it blocks

library(bayesim)

# =============================================================================
# 1. Error Handling - Error classes are properly thrown
# =============================================================================

describe("Error Handling", {
  describe("bayesim_config_error()", {
    it("creates condition with correct class hierarchy", {
      err <- bayesim_config_error("Invalid configuration")
      expect_s3_class(err, "bayesim_config_error")
      expect_s3_class(err, "bayesim_error")
      expect_s3_class(err, "error")
      expect_s3_class(err, "condition")
    })

    it("is properly classified as a fatal error", {
      err <- bayesim_config_error("Invalid configuration")
      expect_true(is_fatal_error(err))
      expect_false(is_recoverable_error(err))
    })

    it("works with tryCatch using bayesim_error handler", {
      caught <- FALSE
      tryCatch(
        stop(bayesim_config_error("Test error")),
        bayesim_error = function(e) {
          caught <<- TRUE
          expect_s3_class(e, "bayesim_config_error")
        }
      )
      expect_true(caught)
    })

    it("works with tryCatch using bayesim_config_error handler", {
      caught <- FALSE
      tryCatch(
        stop(bayesim_config_error("Test error")),
        bayesim_config_error = function(e) {
          caught <<- TRUE
        }
      )
      expect_true(caught)
    })
  })

  describe("bayesim_checkpoint_error()", {
    it("creates condition with correct class hierarchy", {
      err <- bayesim_checkpoint_error("Checkpoint write failed")
      expect_s3_class(err, "bayesim_checkpoint_error")
      expect_s3_class(err, "bayesim_error")
      expect_s3_class(err, "error")
      expect_s3_class(err, "condition")
    })

    it("is properly classified as a fatal error", {
      err <- bayesim_checkpoint_error("Checkpoint failed")
      expect_true(is_fatal_error(err))
      expect_false(is_recoverable_error(err))
    })

    it("works with tryCatch", {
      result <- tryCatch(
        stop(bayesim_checkpoint_error("Checkpoint failed")),
        bayesim_checkpoint_error = function(e) {
          "caught_checkpoint_error"
        },
        error = function(e) {
          "caught_generic_error"
        }
      )
      expect_equal(result, "caught_checkpoint_error")
    })
  })

  describe("write_json_atomic()", {
    it("throws bayesim_checkpoint_error when file.rename fails", {
      tmp_file <- tempfile(fileext = ".json")
      on.exit(if (file.exists(tmp_file)) file.remove(tmp_file))

      # Mock file.rename to return FALSE
      local_mocked_bindings(
        file.rename = function(from, to) FALSE,
        .package = "base"
      )

      expect_error(
        write_json_atomic(list(x = 1), tmp_file),
        class = "bayesim_checkpoint_error"
      )
    })

    it("throws error that is catchable as bayesim_checkpoint_error", {
      tmp_file <- tempfile(fileext = ".json")
      on.exit(if (file.exists(tmp_file)) file.remove(tmp_file))

      # Mock file.rename to return FALSE
      local_mocked_bindings(
        file.rename = function(from, to) FALSE,
        .package = "base"
      )

      caught_class <- NULL
      tryCatch(
        write_json_atomic(list(x = 1), tmp_file),
        bayesim_checkpoint_error = function(e) {
          caught_class <<- class(e)
        }
      )

      expect_true("bayesim_checkpoint_error" %in% caught_class)
    })

    it("succeeds when file.rename returns TRUE", {
      tmp_file <- tempfile(fileext = ".json")
      on.exit(if (file.exists(tmp_file)) file.remove(tmp_file))

      # Should succeed without mocking
      expect_no_error(
        write_json_atomic(list(x = 1, y = "test"), tmp_file)
      )
      expect_true(file.exists(tmp_file))
    })
  })
})

# =============================================================================
# 2. Result Constructor Defaults
# =============================================================================

describe("Result Constructor Defaults", {
  describe("new_fit_result()", {
    it("defaults success to TRUE", {
      result <- new_fit_result()
      expect_true(result$success)
    })

    it("defaults timing components to 0", {
      result <- new_fit_result()
      expect_equal(result$timing$total, 0)
      expect_equal(result$timing$warmup, 0)
      expect_equal(result$timing$sample, 0)
    })

    it("defaults warnings to empty character vector", {
      result <- new_fit_result()
      expect_equal(result$warnings, character())
      expect_type(result$warnings, "character")
    })

    it("defaults diagnostics to empty list", {
      result <- new_fit_result()
      expect_equal(result$diagnostics, list())
      expect_type(result$diagnostics, "list")
    })

    it("defaults fit to NULL", {
      result <- new_fit_result()
      expect_null(result$fit)
    })

    it("defaults draws to NULL", {
      result <- new_fit_result()
      expect_null(result$draws)
    })

    it("defaults error to NULL when success is TRUE", {
      result <- new_fit_result()
      expect_null(result$error)
    })

    it("handles explicit NULL timing by defaulting all components to 0", {
      result <- new_fit_result(timing = NULL)
      expect_equal(result$timing$total, 0)
      expect_equal(result$timing$warmup, 0)
      expect_equal(result$timing$sample, 0)
    })

    it("handles partial timing list", {
      result <- new_fit_result(timing = list(total = 5.0))
      expect_equal(result$timing$total, 5.0)
      expect_equal(result$timing$warmup, 0)
      expect_equal(result$timing$sample, 0)
    })
  })

  describe("new_task_result()", {
    it("defaults status to 'success'", {
      result <- new_task_result(task_id = "test_task", metrics = list(x = 1))
      expect_equal(result$status, "success")
    })

    it("defaults timing$total to 0", {
      result <- new_task_result(task_id = "test_task", metrics = list(x = 1))
      expect_equal(result$timing$total, 0)
    })

    it("defaults warnings to empty character vector", {
      result <- new_task_result(task_id = "test_task", metrics = list(x = 1))
      expect_equal(result$warnings, character())
    })

    it("handles NULL timing by defaulting total to 0", {
      result <- new_task_result(
        task_id = "test_task",
        metrics = list(x = 1),
        timing = NULL
      )
      expect_equal(result$timing$total, 0)
    })

    it("requires metrics when status is 'success'", {
      # This should work - default status is 'success' with metrics
      result <- new_task_result(
        task_id = "test_task",
        metrics = list(rmse = 0.1)
      )
      expect_equal(result$status, "success")
      expect_equal(result$metrics$rmse, 0.1)
    })

    it("allows failed status with error but no metrics", {
      result <- new_task_result(
        task_id = "test_task",
        status = "failed",
        error = list(error_class = "test", error_message = "failed")
      )
      expect_equal(result$status, "failed")
      expect_null(result$metrics)
    })
  })

  describe("new_simulation_result()", {
    it("defaults timing$total to 0", {
      result <- new_simulation_result()
      expect_equal(result$timing$total, 0)
    })

    it("defaults task_results to empty list", {
      result <- new_simulation_result()
      expect_equal(result$task_results, list())
    })

    it("defaults summary to empty data.frame", {
      result <- new_simulation_result()
      expect_true(is.data.frame(result$summary))
      expect_equal(nrow(result$summary), 0)
    })

    it("defaults errors to empty data.frame", {
      result <- new_simulation_result()
      expect_true(is.data.frame(result$errors))
      expect_equal(nrow(result$errors), 0)
    })

    it("defaults checkpoint_path to NULL", {
      result <- new_simulation_result()
      expect_null(result$checkpoint_path)
    })

    it("handles NULL timing by defaulting total to 0", {
      result <- new_simulation_result(timing = NULL)
      expect_equal(result$timing$total, 0)
    })
  })

  describe("Empty/zero-length defaults validation", {
    it("validates fit_result with empty warnings and diagnostics", {
      result <- new_fit_result(
        warnings = character(),
        diagnostics = list()
      )
      expect_no_error(validate_bayesim_fit_result(result))
    })

    it("validates task_result with empty warnings", {
      result <- new_task_result(
        task_id = "test",
        metrics = list(x = 1),
        warnings = character()
      )
      expect_no_error(validate_bayesim_task_result(result))
    })

    it("validates simulation_result with empty collections", {
      result <- new_simulation_result(
        task_results = list(),
        summary = data.frame(),
        errors = data.frame()
      )
      expect_no_error(validate_bayesim_simulation_result(result))
    })
  })
})

# =============================================================================
# 3. MockFitter Contract Compliance
# =============================================================================

describe("MockFitter Contract Compliance", {
  fitter <- MockFitter()
  data_bundle <- list(
    train = data.frame(y = 1:10, x = 1:10),
    test = NULL,
    response = "y",
    true_params = c(beta = 1.0),
    vars_of_interest = "beta",
    references = c(beta = 0),
    meta = list()
  )
  fit_spec <- list(model = "baseline")
  task_ctx <- list(
    task_id = "d001_f001_r00001",
    data_idx = 1,
    fit_idx = 1,
    rep_idx = 1
  )

  describe("fit() return value", {
    it("returns bayesim_fit_result (not just classed list)", {
      result <- fit(
        fitter,
        data_bundle,
        fit_spec,
        seed = 123L,
        task_ctx = task_ctx
      )
      expect_true(is_bayesim_fit_result(result))
      expect_s3_class(result, "bayesim_fit_result")
    })

    it("returns result with all required fields", {
      result <- fit(
        fitter,
        data_bundle,
        fit_spec,
        seed = 123L,
        task_ctx = task_ctx
      )

      # Check all required fields exist
      expect_true("success" %in% names(result))
      expect_true("fit" %in% names(result))
      expect_true("draws" %in% names(result))
      expect_true("diagnostics" %in% names(result))
      expect_true("timing" %in% names(result))
      expect_true("warnings" %in% names(result))
      expect_true("error" %in% names(result))
    })

    it("returns success=TRUE for valid input", {
      result <- fit(
        fitter,
        data_bundle,
        fit_spec,
        seed = 123L,
        task_ctx = task_ctx
      )
      expect_true(result$success)
    })

    it("returns draws matrix with colnames", {
      result <- fit(
        fitter,
        data_bundle,
        fit_spec,
        seed = 123L,
        task_ctx = task_ctx
      )

      expect_true(is.matrix(result$draws))
      expect_false(is.null(colnames(result$draws)))
      expect_true(length(colnames(result$draws)) > 0)
    })

    it("returns timing with total, warmup, sample components", {
      result <- fit(
        fitter,
        data_bundle,
        fit_spec,
        seed = 123L,
        task_ctx = task_ctx
      )

      expect_true(is.list(result$timing))
      expect_true("total" %in% names(result$timing))
      expect_true("warmup" %in% names(result$timing))
      expect_true("sample" %in% names(result$timing))

      expect_type(result$timing$total, "double")
      expect_type(result$timing$warmup, "double")
      expect_type(result$timing$sample, "double")
    })

    it("returns diagnostics as named list", {
      result <- fit(
        fitter,
        data_bundle,
        fit_spec,
        seed = 123L,
        task_ctx = task_ctx
      )

      expect_true(is.list(result$diagnostics))
      expect_true(length(result$diagnostics) > 0)
      expect_true(!is.null(names(result$diagnostics)))
    })

    it("returns warnings as character vector", {
      result <- fit(
        fitter,
        data_bundle,
        fit_spec,
        seed = 123L,
        task_ctx = task_ctx
      )

      expect_type(result$warnings, "character")
    })

    it("returns NULL error when successful", {
      result <- fit(
        fitter,
        data_bundle,
        fit_spec,
        seed = 123L,
        task_ctx = task_ctx
      )
      expect_null(result$error)
    })
  })

  describe("draws matrix structure", {
    it("has colnames for all parameters", {
      result <- fit(
        fitter,
        data_bundle,
        fit_spec,
        seed = 123L,
        task_ctx = task_ctx
      )

      param_names <- colnames(result$draws)
      expect_equal(length(param_names), ncol(result$draws))
      expect_false(anyNA(param_names))
      expect_false(any(param_names == ""))
    })
  })

  describe("validation", {
    it("passes validate_bayesim_fit_result", {
      result <- fit(
        fitter,
        data_bundle,
        fit_spec,
        seed = 123L,
        task_ctx = task_ctx
      )
      expect_no_error(validate_bayesim_fit_result(result))
    })
  })
})

# =============================================================================
# 4. Serialization Safety
# =============================================================================

describe("Serialization Safety", {
  describe("capture_fitter_spec()", {
    it("returns NA for NULL input", {
      result <- capture_fitter_spec(NULL)
      expect_true(is.na(result))
    })

    it("never digests S7 objects - returns plain list", {
      fitter <- MockFitter()
      result <- capture_fitter_spec(fitter)

      # Result should be a plain list, not an S7 object
      expect_type(result, "list")
      expect_false(S7::S7_inherits(result))
    })

    it("returns plain list with class, name for S7 fitter", {
      fitter <- MockFitter()
      result <- capture_fitter_spec(fitter)

      expect_true("class" %in% names(result))
      expect_true("properties" %in% names(result))
      expect_true(grepl("MockFitter", result$class, fixed = TRUE))
      expect_equal(result$properties$name, "mock")
    })

    it("includes supports_* properties when available", {
      fitter <- MockFitter()
      result <- capture_fitter_spec(fitter)

      expect_true("supports_predictions" %in% names(result$properties))
      expect_true("supports_log_lik" %in% names(result$properties))
      expect_true("supports_loo" %in% names(result$properties))
    })

    it("returns NA for non-S7 objects", {
      result <- capture_fitter_spec(list(name = "test"))
      expect_true(is.na(result))
    })
  })

  describe("capture_metrics_spec()", {
    it("returns empty list for NULL input", {
      result <- capture_metrics_spec(NULL)
      expect_equal(result, list())
    })

    it("never digests S7 objects - returns list of plain lists", {
      metrics <- list(rmse_metric())
      result <- capture_metrics_spec(metrics)

      # Result should be a list of lists
      expect_type(result, "list")
      expect_false(S7::S7_inherits(result))

      # Each element should also be plain
      for (elem in result) {
        expect_type(elem, "list")
        expect_true("properties" %in% names(elem))
      }
    })

    it("captures S7 metric specs as plain metadata", {
      metrics <- list(rmse_metric(), bias_metric())
      result <- capture_metrics_spec(metrics)

      expect_type(result, "list")
      expect_equal(length(result), 2)
      expect_equal(result[[1]]$properties$name, "rmse")
      expect_equal(result[[2]]$properties$name, "bias")
    })

    it("handles unknown types gracefully with NA", {
      metrics <- list("unknown_type")
      result <- capture_metrics_spec(metrics)

      expect_type(result, "list")
      expect_equal(length(result), 1)
      expect_true(is.na(result[[1]]))
    })
  })
})

# =============================================================================
# 5. Consistent Separator in flatten_with_prefix
# =============================================================================

describe("flatten_with_prefix()", {
  it("uses __ for all separators", {
    x <- list(a = c(x = 1, y = 2), b = c(p = 3, q = 4))
    result <- flatten_with_prefix(x, "param")

    # All flattened names should use __ separator
    flat_names <- grep("^param__", names(result), value = TRUE)
    expect_true(length(flat_names) > 0)

    # Check the pattern is prefix__name__subname
    expect_true("param__a__x" %in% names(result))
    expect_true("param__a__y" %in% names(result))
    expect_true("param__b__p" %in% names(result))
    expect_true("param__b__q" %in% names(result))
  })

  it("does not use single underscore between name and subname", {
    x <- list(a = c(x = 1, y = 2))
    result <- flatten_with_prefix(x, "param")

    # Should NOT have param__a_x (single underscore)
    expect_false("param__a_x" %in% names(result))
    expect_false("param__a_y" %in% names(result))

    # Should have param__a__x (double underscore)
    expect_true("param__a__x" %in% names(result))
    expect_true("param__a__y" %in% names(result))
  })

  it("preserves scalar values without modification", {
    x <- list(scalar_val = 42, nested = c(a = 1, b = 2))
    result <- flatten_with_prefix(x, "param")

    expect_true("scalar_val" %in% names(result))
    expect_equal(result$scalar_val, 42)
  })

  it("preserves unnamed numeric vectors (length > 1) as-is", {
    x <- list(unnamed = c(1, 2, 3))
    result <- flatten_with_prefix(x, "param")

    # Unnamed vectors should not be flattened
    expect_true("unnamed" %in% names(result))
    expect_equal(result$unnamed, c(1, 2, 3))
  })

  it("handles single element named vectors as scalars", {
    x <- list(single = c(a = 1))
    result <- flatten_with_prefix(x, "param")

    # Single element named vectors should be treated as scalars (not flattened)
    expect_true("single" %in% names(result))
    expect_equal(result$single, c(a = 1))
  })

  it("handles empty list gracefully", {
    result <- flatten_with_prefix(list(), "param")
    expect_type(result, "list")
    expect_equal(length(result), 0)
  })
})

# =============================================================================
# 6. NULL Safety in MockFitter
# =============================================================================

describe("NULL Safety in MockFitter", {
  fitter <- MockFitter()

  describe("fit() with NULL data_bundle", {
    it("throws error for NULL data_bundle$train", {
      expect_error(
        fit(fitter, list(train = NULL), list(), seed = 123L, task_ctx = list()),
        class = "bayesim_contract_error"
      )
    })
  })

  describe("extract_draws() with NULL components", {
    it("handles NULL seed in fit_result by using default", {
      # Create a mock fit object with NULL seed
      mock_fit <- list(
        fitter = "mock",
        data_bundle = list(train = data.frame(y = 1:10, x = 1:10)),
        n_obs = 10,
        seed = NULL
      )

      # Should not error, uses default seed 12345L
      result <- extract_draws(fitter, mock_fit)
      expect_true(is.matrix(result))
      expect_false(is.null(colnames(result)))
    })

    it("handles NULL fit_result by extracting from list structure", {
      # When called internally from fit(), it passes a list not a bayesim_fit_result
      fit_obj <- list(
        seed = 42L,
        data_bundle = list(train = data.frame(y = 1:5)),
        n_obs = 5
      )

      result <- extract_draws(fitter, fit_obj)
      expect_true(is.matrix(result))
      expect_true(ncol(result) > 0)
    })
  })

  describe("predict_fit() with NULL newdata", {
    it("handles NULL newdata by using training data", {
      data_bundle <- list(
        train = data.frame(y = 1:10, x = 1:10),
        test = NULL,
        response = "y"
      )

      fit_result <- fit(
        fitter,
        data_bundle,
        list(model = "baseline"),
        seed = 123L,
        task_ctx = list(task_id = "test")
      )

      # Should not error when newdata is NULL
      predictions <- predict_fit(fitter, fit_result, newdata = NULL)
      expect_type(predictions, "list")
      expect_true("predicted_mean" %in% names(predictions))
    })
  })

  describe("log_lik() with NULL newdata", {
    it("handles NULL newdata by using training data", {
      data_bundle <- list(
        train = data.frame(y = 1:10, x = 1:10),
        test = NULL,
        response = "y"
      )

      fit_result <- fit(
        fitter,
        data_bundle,
        list(model = "baseline"),
        seed = 123L,
        task_ctx = list(task_id = "test")
      )

      # Should not error when newdata is NULL
      ll <- log_lik(fitter, fit_result, newdata = NULL)
      expect_true(is.matrix(ll))
      expect_equal(nrow(ll), 10) # n_obs from training data
    })
  })

  describe("diagnostics() with NULL/empty diagnostics", {
    it("returns default diagnostics for fit_result with empty diagnostics", {
      data_bundle <- list(
        train = data.frame(y = 1:10, x = 1:10)
      )

      # Create a fit_result with empty diagnostics
      fit_result <- list(
        diagnostics = list(),
        fit = list(seed = 123L, n_obs = 10)
      )
      class(fit_result) <- "bayesim_fit_result"

      result <- diagnostics(fitter, fit_result)
      expect_type(result, "list")
      expect_true("rhat_max" %in% names(result))
    })
  })
})
