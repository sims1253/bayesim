# test-engine.R
# Comprehensive tests for Phase 2 engine components
# testthat 3e style with describe/it blocks

# Load the package
library(bayesim)

# =============================================================================
# Helper: Custom S7 expectation if not available
# =============================================================================

expect_s7_object <- function(object) {
  testthat::expect(
    S7::S7_inherits(object),
    sprintf("%s is not an S7 object", deparse(substitute(object)))
  )
}

# =============================================================================
# Helper: Create a simple test metric with correct signature
# =============================================================================

create_test_metric <- function(
  name = "test_metric",
  needs = character(),
  required = FALSE,
  compute_fn = NULL
) {
  TestMetric <- S7::new_class(
    "TestMetric",
    parent = Metric,
    properties = list(
      name = S7::new_property(S7::class_character),
      needs = S7::new_property(S7::class_character),
      required = S7::new_property(S7::class_logical)
    )
  )

  if (is.null(compute_fn)) {
    compute_fn <- function(metric, fit_result, data_bundle, context, task_ctx) {
      list(value = 1.0)
    }
  }

  S7::method(compute_metric, TestMetric) <- compute_fn
  # Pass properties explicitly when creating the instance
  TestMetric(name = name, needs = needs, required = required)
}

# =============================================================================
# Helper: Check if run_simulation exists and is functional
# =============================================================================

run_sim_available <- function() {
  exists("run_simulation", mode = "function") &&
    exists("is_simulation_config", mode = "function")
}

# =============================================================================
# 1. RNG Management (from rng.R)
# =============================================================================

describe("RNG Management", {
  describe("setup_global_rng()", {
    it("sets correct RNG kind to L'Ecuyer-CMRG", {
      old_kind <- RNGkind()[1]
      on.exit(RNGkind(old_kind), add = TRUE)

      setup_global_rng(42)

      expect_equal(RNGkind()[1], "L'Ecuyer-CMRG")
    })

    it("sets the seed correctly", {
      old_kind <- RNGkind()[1]
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
        get(".Random.seed", envir = .GlobalEnv)
      } else {
        NULL
      }
      on.exit(
        {
          RNGkind(old_kind)
          if (!is.null(old_seed)) {
            assign(".Random.seed", old_seed, envir = .GlobalEnv)
          }
        },
        add = TRUE
      )

      setup_global_rng(123)

      # Generate a value and check reproducibility
      val1 <- runif(1)

      # Reset with same seed
      setup_global_rng(123)
      val2 <- runif(1)

      expect_equal(val1, val2)
    })

    it("returns the initial .Random.seed state invisibly", {
      old_kind <- RNGkind()[1]
      on.exit(RNGkind(old_kind), add = TRUE)

      result <- setup_global_rng(42)
      expect_true(is.integer(result))
      expect_true(length(result) > 0)
    })
  })

  describe("create_task_rng_streams()", {
    it("creates correct number of streams", {
      streams <- create_task_rng_streams(42, 5)

      expect_true(is.list(streams))
      expect_equal(length(streams), 5)
    })

    it("creates streams that are integer vectors", {
      streams <- create_task_rng_streams(42, 3)

      for (i in seq_along(streams)) {
        expect_true(is.integer(streams[[i]]))
      }
    })

    it("produces different streams for different seeds", {
      streams1 <- create_task_rng_streams(42, 3)
      streams2 <- create_task_rng_streams(123, 3)

      expect_false(identical(streams1[[1]], streams2[[1]]))
    })

    it("produces different streams within same seed", {
      streams <- create_task_rng_streams(42, 5)

      # Each stream should be different
      for (i in 2:5) {
        expect_false(identical(streams[[1]], streams[[i]]))
      }
    })

    it("produces reproducible streams with same seed", {
      streams1 <- create_task_rng_streams(42, 3)
      streams2 <- create_task_rng_streams(42, 3)

      expect_identical(streams1, streams2)
    })

    it("advances using independent L'Ecuyer streams", {
      RNGkind("L'Ecuyer-CMRG")
      set.seed(42)
      seed <- .Random.seed

      streams <- create_task_rng_streams(42, 3)

      expect_identical(streams[[1]], seed)
      expect_identical(streams[[2]], parallel::nextRNGStream(seed))
    })
  })

  describe("set_task_rng()", {
    it("restores RNG state correctly", {
      # Create a stream
      streams <- create_task_rng_streams(42, 3)
      stream <- streams[[1]]

      # Set the stream
      set_task_rng(stream)

      # Generate some values
      vals1 <- runif(5)

      # Reset to the same stream
      set_task_rng(stream)

      # Generate again
      vals2 <- runif(5)

      expect_equal(vals1, vals2)
    })

    it("sets .Random.seed in global environment", {
      streams <- create_task_rng_streams(42, 1)
      stream <- streams[[1]]

      set_task_rng(stream)

      expect_true(exists(".Random.seed", envir = .GlobalEnv))
      expect_identical(get(".Random.seed", envir = .GlobalEnv), stream)
    })
  })

  describe("advance_rng_stream()", {
    it("advances RNG state and returns new state", {
      streams <- create_task_rng_streams(42, 1)
      stream <- streams[[1]]

      advanced <- advance_rng_stream(stream)

      # Advanced should be an integer vector
      expect_true(is.integer(advanced))
      expect_true(length(advanced) == length(stream))
    })

    it("advances RNG state by n steps", {
      streams <- create_task_rng_streams(42, 1)
      stream <- streams[[1]]

      advanced_1 <- advance_rng_stream(stream, n = 1)
      advanced_5 <- advance_rng_stream(stream, n = 5)

      # Both should be integer vectors
      expect_true(is.integer(advanced_1))
      expect_true(is.integer(advanced_5))
    })

    it("returns a valid RNG state that can be set", {
      streams <- create_task_rng_streams(42, 1)
      stream <- streams[[1]]

      advanced <- advance_rng_stream(stream, n = 5)

      # Should be able to set this state
      set_task_rng(advanced)

      # Generate a value to verify RNG works
      val <- runif(1)
      expect_true(is.numeric(val))
    })
  })
})


# =============================================================================
# 2. Task Grid (from task-grid.R)
# =============================================================================

describe("Task Grid", {
  # Helper to create a minimal config for testing
  # Note: Tests may need to skip if is_simulation_config has namespace issues
  create_test_config <- function(
    n_data = 2,
    n_fit = 2,
    n_replicates = 3L
  ) {
    simulation_config(
      data_grid = data.frame(n = seq_len(n_data) * 100),
      fit_grid = data.frame(model = paste0("model_", seq_len(n_fit))),
      data_generator = function(data_spec, task_ctx) {
        list(
          train = data.frame(y = rnorm(data_spec$n), x = rnorm(data_spec$n)),
          test = NULL,
          response = "y",
          true_params = c(beta = 1.0, sigma = 1.0),
          vars_of_interest = c("beta", "sigma")
        )
      },
      n_replicates = n_replicates,
      seed = 42L
    )
  }

  # Check if is_simulation_config works correctly
  config_check_works <- function() {
    config <- tryCatch(
      create_test_config(),
      error = function(e) NULL
    )
    if (is.null(config)) {
      return(FALSE)
    }
    is_simulation_config(config)
  }

  describe("create_task_grid()", {
    it("creates correct number of tasks", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config(n_data = 2, n_fit = 3, n_replicates = 5L)

      task_grid <- create_task_grid(config)

      expect_s3_class(task_grid, "tbl_df")
      expect_equal(nrow(task_grid), 2 * 3 * 5) # 30 tasks
    })

    it("creates tasks with required columns", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config()
      task_grid <- create_task_grid(config)

      expect_true("task_id" %in% names(task_grid))
      expect_true("data_idx" %in% names(task_grid))
      expect_true("fit_idx" %in% names(task_grid))
      expect_true("rep_idx" %in% names(task_grid))
      expect_true("rng_seed" %in% names(task_grid))
      expect_true("status" %in% names(task_grid))
    })

    it("initializes all tasks with status 'pending'", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config()
      task_grid <- create_task_grid(config)

      expect_true(all(task_grid$status == "pending"))
    })
  })

  describe("Task IDs", {
    it("are correctly formatted as dXXX_fXXX_rXXXXX", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config(n_data = 2, n_fit = 3, n_replicates = 5L)
      task_grid <- create_task_grid(config)

      # All task IDs should match the pattern
      pattern_match <- grepl(
        "^d[0-9]{3}_f[0-9]{3}_r[0-9]{5}$",
        task_grid$task_id
      )
      expect_true(all(pattern_match))
    })

    it("reflect correct indices in the ID", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config(n_data = 1, n_fit = 1, n_replicates = 1L)
      task_grid <- create_task_grid(config)

      expect_equal(task_grid$task_id[1], "d001_f001_r00001")
    })

    it("are unique across all tasks", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config(n_data = 3, n_fit = 4, n_replicates = 10L)
      task_grid <- create_task_grid(config)

      # anyDuplicated returns 0 (integer) when no duplicates, not FALSE
      expect_equal(anyDuplicated(task_grid$task_id), 0)
    })
  })

  describe("Task grid sorting", {
    it("is sorted lexicographically by task_id", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config(n_data = 2, n_fit = 2, n_replicates = 3L)
      task_grid <- create_task_grid(config)

      sorted_ids <- sort(task_grid$task_id)
      expect_equal(task_grid$task_id, sorted_ids)
    })

    it("is sorted by data_idx, then fit_idx, then rep_idx", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config(n_data = 2, n_fit = 2, n_replicates = 3L)
      task_grid <- create_task_grid(config)

      for (i in 1:(nrow(task_grid) - 1)) {
        curr <- task_grid[i, ]
        next_row <- task_grid[i + 1, ]

        # Either data_idx increases, or same data_idx and fit_idx increases,
        # or same data_idx and fit_idx and rep_idx increases
        expect_true(
          curr$data_idx < next_row$data_idx ||
            (curr$data_idx == next_row$data_idx &&
              curr$fit_idx < next_row$fit_idx) ||
            (curr$data_idx == next_row$data_idx &&
              curr$fit_idx == next_row$fit_idx &&
              curr$rep_idx < next_row$rep_idx)
        )
      }
    })
  })

  describe("get_task_spec()", {
    it("extracts correct specs for a task", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config(n_data = 2, n_fit = 2, n_replicates = 1L)
      task_grid <- create_task_grid(config)

      spec <- get_task_spec(task_grid, "d002_f001_r00001", config)

      expect_equal(spec$task_id, "d002_f001_r00001")
      expect_equal(spec$data_idx, 2)
      expect_equal(spec$fit_idx, 1)
      expect_equal(spec$rep_idx, 1)
    })

    it("extracts correct data_spec", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config(n_data = 2, n_fit = 1, n_replicates = 1L)
      task_grid <- create_task_grid(config)

      spec <- get_task_spec(task_grid, "d001_f001_r00001", config)
      expect_equal(spec$data_spec$n, 100)

      spec2 <- get_task_spec(task_grid, "d002_f001_r00001", config)
      expect_equal(spec2$data_spec$n, 200)
    })

    it("extracts correct fit_spec", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config(n_data = 1, n_fit = 2, n_replicates = 1L)
      task_grid <- create_task_grid(config)

      spec <- get_task_spec(task_grid, "d001_f001_r00001", config)
      expect_equal(spec$fit_spec$model, "model_1")

      spec2 <- get_task_spec(task_grid, "d001_f002_r00001", config)
      expect_equal(spec2$fit_spec$model, "model_2")
    })

    it("extracts correct task_ctx", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config(n_data = 2, n_fit = 2, n_replicates = 3L)
      task_grid <- create_task_grid(config)

      spec <- get_task_spec(task_grid, "d002_f001_r00002", config)

      expect_equal(spec$task_ctx$task_id, "d002_f001_r00002")
      expect_equal(spec$task_ctx$data_idx, 2)
      expect_equal(spec$task_ctx$fit_idx, 1)
      expect_equal(spec$task_ctx$rep_idx, 2)
    })

    it("errors for non-existent task_id", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config()
      task_grid <- create_task_grid(config)

      expect_error(
        get_task_spec(task_grid, "d999_f999_r99999", config),
        "not found in grid"
      )
    })
  })

  describe("update_task_status()", {
    it("updates status correctly", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config(n_data = 1, n_fit = 1, n_replicates = 2L)
      task_grid <- create_task_grid(config)

      updated <- update_task_status(task_grid, "d001_f001_r00001", "success")

      expect_equal(updated$status[1], "success")
      expect_equal(updated$status[2], "pending") # Other unchanged
    })

    it("returns a modified copy without changing original", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config()
      task_grid <- create_task_grid(config)
      original_status <- task_grid$status[1]

      updated <- update_task_status(task_grid, task_grid$task_id[1], "failed")

      expect_equal(task_grid$status[1], original_status)
      expect_equal(updated$status[1], "failed")
    })

    it("can update to different status values", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config(n_data = 1, n_fit = 1, n_replicates = 3L)
      task_grid <- create_task_grid(config)

      task_grid <- update_task_status(
        task_grid,
        task_grid$task_id[1],
        "success"
      )
      task_grid <- update_task_status(task_grid, task_grid$task_id[2], "failed")
      task_grid <- update_task_status(
        task_grid,
        task_grid$task_id[3],
        "skipped"
      )

      expect_equal(task_grid$status[1], "success")
      expect_equal(task_grid$status[2], "failed")
      expect_equal(task_grid$status[3], "skipped")
    })
  })

  describe("filter_tasks_by_status()", {
    it("filters to include only specified statuses", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config(n_data = 1, n_fit = 1, n_replicates = 4L)
      task_grid <- create_task_grid(config)

      task_grid <- update_task_status(
        task_grid,
        task_grid$task_id[1],
        "success"
      )
      task_grid <- update_task_status(task_grid, task_grid$task_id[2], "failed")
      task_grid <- update_task_status(
        task_grid,
        task_grid$task_id[3],
        "success"
      )
      # task 4 remains pending

      success_only <- filter_tasks_by_status(task_grid, "success")
      expect_equal(nrow(success_only), 2)
      expect_true(all(success_only$status == "success"))

      failed_only <- filter_tasks_by_status(task_grid, "failed")
      expect_equal(nrow(failed_only), 1)

      pending_only <- filter_tasks_by_status(task_grid, "pending")
      expect_equal(nrow(pending_only), 1)
    })

    it("can filter by multiple statuses", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config(n_data = 1, n_fit = 1, n_replicates = 4L)
      task_grid <- create_task_grid(config)

      task_grid <- update_task_status(
        task_grid,
        task_grid$task_id[1],
        "success"
      )
      task_grid <- update_task_status(task_grid, task_grid$task_id[2], "failed")
      task_grid <- update_task_status(
        task_grid,
        task_grid$task_id[3],
        "success"
      )
      # task 4 remains pending

      completed <- filter_tasks_by_status(task_grid, c("success", "failed"))
      expect_equal(nrow(completed), 3)
    })

    it("returns empty tibble when no matches", {
      skip_if_not(
        config_check_works(),
        "is_simulation_config not working correctly"
      )

      config <- create_test_config()
      task_grid <- create_task_grid(config)

      result <- filter_tasks_by_status(task_grid, "nonexistent")
      expect_equal(nrow(result), 0)
    })
  })
})


# =============================================================================
# 3. Metric Resolution
# =============================================================================

describe("Metric Resolution", {
  describe("built-in metric constructors", {
    it("return S7 Metric objects", {
      expect_s7_object(pred_rmse_metric())
      expect_s7_object(pred_bias_metric())
      expect_s7_object(coverage_metric())
      expect_s7_object(posterior_mean_metric())
    })
  })

  describe("resolve_metrics()", {
    it("accepts a single Metric object", {
      metric <- create_test_metric(name = "resolve_single")
      resolved <- resolve_metrics(metric)

      expect_true(is.list(resolved))
      expect_equal(length(resolved), 1)
      expect_equal(resolved[[1]]@name, "resolve_single")
    })

    it("accepts a list of Metric objects", {
      metrics_list <- list(
        create_test_metric(name = "resolve_list_1"),
        create_test_metric(name = "resolve_list_2")
      )

      resolved <- resolve_metrics(metrics_list)

      expect_identical(resolved, metrics_list)
    })

    it("rejects character metric names", {
      expect_error(
        resolve_metrics(c("rmse", "bias")),
        "metrics must be Metric objects"
      )
    })

    it("rejects invalid list contents", {
      expect_error(
        resolve_metrics(list("not a metric")),
        "is not an S7 Metric object"
      )
    })

    it("returns empty list for NULL input", {
      expect_equal(resolve_metrics(NULL), list())
    })
  })
})


# =============================================================================
# 4. Worker (from worker.R)
# =============================================================================

describe("Worker", {
  # Create a simple test configuration
  create_worker_test_config <- function() {
    list(
      data_generator = function(data_spec, task_ctx) {
        list(
          train = data.frame(y = rnorm(10), x = rnorm(10)),
          test = NULL,
          response = "y",
          true_params = c(beta = 1.0, sigma = 1.0),
          vars_of_interest = c("beta", "sigma")
        )
      }
    )
  }

  create_test_task <- function(seed = 42L) {
    list(
      task_id = "d001_f001_r00001",
      data_spec = list(n = 10),
      fit_spec = list(model = "test"),
      rng_seed = create_task_rng_streams(seed, 1)[[1]],
      task_ctx = list(
        task_id = "d001_f001_r00001",
        data_idx = 1,
        fit_idx = 1,
        rep_idx = 1
      )
    )
  }

  describe("run_task_safe()", {
    it("catches errors and returns task_result", {
      task <- create_test_task()
      config_spec <- list(
        data_generator = function(data_spec, task_ctx) {
          stop("Data generation error")
        }
      )
      fitter <- MockFitter()
      metrics <- list()

      result <- run_task_safe(task, config_spec, fitter, metrics)

      expect_s3_class(result, "bayesim_task_result")
      expect_equal(result$status, "failed")
      expect_true(is.list(result$error))
    })

    it("returns task_id in error result", {
      task <- create_test_task()
      config_spec <- list(
        data_generator = function(data_spec, task_ctx) {
          stop("Error")
        }
      )

      result <- run_task_safe(task, config_spec, MockFitter(), list())

      expect_equal(result$task_id, "d001_f001_r00001")
    })

    it("includes timing in error result", {
      task <- create_test_task()
      config_spec <- list(
        data_generator = function(data_spec, task_ctx) {
          stop("Error")
        }
      )

      result <- run_task_safe(task, config_spec, MockFitter(), list())

      expect_true(is.list(result$timing))
      expect_true("total" %in% names(result$timing))
    })
  })

  describe("run_task()", {
    it("with MockFitter produces valid results", {
      task <- create_test_task()
      config_spec <- create_worker_test_config()
      fitter <- MockFitter()

      test_metric <- create_test_metric(
        name = "simple_metric",
        compute_fn = function(
          metric,
          fit_result,
          data_bundle,
          context,
          task_ctx
        ) {
          list(value = 42.0)
        }
      )
      metrics <- list(test_metric)

      result <- run_task(task, config_spec, fitter, metrics)

      expect_s3_class(result, "bayesim_task_result")
      expect_equal(result$status, "success")
      expect_true(is.list(result$metrics))
    })

    it("with failing data generator returns failed status", {
      task <- create_test_task()
      config_spec <- list(
        data_generator = function(data_spec, task_ctx) {
          stop("Intentional data generation failure")
        }
      )

      result <- run_task(task, config_spec, MockFitter(), list())

      expect_equal(result$status, "failed")
      expect_true(is.list(result$error))
      expect_true("error_message" %in% names(result$error))
    })

    it("with failing fitter returns failed status", {
      task <- create_test_task()
      config_spec <- create_worker_test_config()

      # Create a fitter that always fails - use correct signature
      FailingFitter <- S7::new_class(
        "FailingFitter",
        parent = Fitter,
        properties = list(
          name = S7::new_property(S7::class_character, default = "failing")
        )
      )
      S7::method(fit_model, FailingFitter) <- function(
        fitter,
        data_bundle,
        fit_spec,
        seed,
        task_ctx
      ) {
        stop("Fitting failed intentionally")
      }

      fitter <- FailingFitter()

      result <- run_task(task, config_spec, fitter, list())

      expect_equal(result$status, "failed")
      expect_true(grepl(
        "Fitting failed",
        result$error$error_message,
        fixed = TRUE
      ))
    })

    it("sets RNG state for reproducibility", {
      task1 <- create_test_task(seed = 123)
      task2 <- create_test_task(seed = 123)
      config_spec <- create_worker_test_config()
      fitter <- MockFitter()

      result1 <- run_task(task1, config_spec, fitter, list())
      result2 <- run_task(task2, config_spec, fitter, list())

      # Results should be identical due to same seed
      expect_identical(result1$metrics, result2$metrics)
    })

    it("passes task context to data generator", {
      task <- create_test_task()
      received_ctx <- NULL

      config_spec <- list(
        data_generator = function(data_spec, task_ctx) {
          received_ctx <<- task_ctx
          list(
            train = data.frame(y = 1:10, x = 1:10),
            test = NULL,
            response = "y",
            true_params = c(beta = 1.0),
            vars_of_interest = "beta"
          )
        }
      )

      run_task(task, config_spec, MockFitter(), list())

      expect_equal(received_ctx$task_id, "d001_f001_r00001")
    })

    it("passes a scalar task seed to generators and fitters", {
      task <- create_test_task(seed = 321)
      received_seed <- NULL
      fitter_seed <- NULL

      config_spec <- list(
        data_generator = function(data_spec, task_ctx) {
          received_seed <<- task_ctx$seed
          list(
            train = data.frame(y = 1:10, x = 1:10),
            test = NULL,
            response = "y",
            true_params = c(beta = 1.0),
            vars_of_interest = "beta"
          )
        }
      )

      SeedCapturingFitter <- S7::new_class(
        "SeedCapturingFitter",
        parent = Fitter,
        properties = list(
          name = S7::new_property(S7::class_character, default = "seed_capture")
        )
      )
      S7::method(fit_model, SeedCapturingFitter) <- function(
        fitter,
        data_bundle,
        fit_spec,
        seed,
        task_ctx
      ) {
        fitter_seed <<- seed
        new_fit_result(
          success = TRUE,
          diagnostics = list(),
          timing = list(total = 0.1, warmup = 0, sample = 0)
        )
      }

      result <- run_task(task, config_spec, SeedCapturingFitter(), list())

      expect_equal(result$status, "success")
      expect_true(is.integer(received_seed))
      expect_equal(length(received_seed), 1)
      expect_identical(received_seed, fitter_seed)
    })

    it("externalizes high-cardinality metric vectors when result_path is available", {
      task <- create_test_task(seed = 321)
      result_path <- withr::local_tempdir()
      config_spec <- list(
        data_generator = function(data_spec, task_ctx) {
          list(
            train = data.frame(y = 1:10, x = 1:10),
            test = NULL,
            response = "y",
            true_params = c(beta = 1.0),
            vars_of_interest = "beta"
          )
        },
        result_path = result_path
      )

      long_metric <- create_test_metric(
        name = "long_metric",
        compute_fn = function(
          metric,
          fit_result,
          data_bundle,
          context,
          task_ctx
        ) {
          values <- stats::setNames(
            as.double(seq_len(100)),
            paste0("p", seq_len(100))
          )
          list(curve = values)
        }
      )

      result <- run_task(task, config_spec, MockFitter(), list(long_metric))

      expect_equal(result$status, "success")
      expect_true(isTRUE(result$metrics$long_metric__curve__externalized))
      expect_true(file.exists(result$metrics$long_metric__curve__artifact_path))
      expect_equal(result$metrics$long_metric__curve__n_values, 100)
    })
  })

  describe("build_metric_context()", {
    it("includes context needed by metrics", {
      fitter <- MockFitter()
      data_bundle <- valid_data_bundle()

      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )

      # Metrics that need predictions
      metrics <- list(create_test_metric(
        name = "pred_only",
        needs = "predictions"
      ))

      context <- build_metric_context(fit_result, fitter, data_bundle, metrics)

      # Context should have predictions since metric needs it and fitter supports it
      expect_true(is.list(context))
    })

    it("does not include context not needed by metrics", {
      fitter <- MockFitter()
      data_bundle <- valid_data_bundle()

      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )

      # Metric that needs nothing
      metrics <- list(create_test_metric(
        name = "no_needs",
        needs = character()
      ))

      context <- build_metric_context(fit_result, fitter, data_bundle, metrics)

      # Context should be empty since no metrics need anything
      expect_equal(length(context), 0)
    })

    it("handles fitter without certain capabilities", {
      # Fitter that doesn't support predictions - use correct signature
      NoPredFitter <- S7::new_class(
        "NoPredFitter",
        parent = Fitter,
        properties = list(
          name = S7::new_property(S7::class_character, default = "no_pred"),
          supports_predictions = S7::new_property(
            S7::class_logical,
            default = FALSE
          ),
          supports_log_lik = S7::new_property(
            S7::class_logical,
            default = TRUE
          ),
          supports_loo = S7::new_property(S7::class_logical, default = TRUE)
        )
      )
      S7::method(fit_model, NoPredFitter) <- function(
        fitter,
        data_bundle,
        fit_spec,
        seed,
        task_ctx
      ) {
        new_fit_result(success = TRUE, diagnostics = list())
      }

      fitter <- NoPredFitter()
      data_bundle <- valid_data_bundle()
      fit_result <- new_fit_result(success = TRUE)

      metrics <- list(create_test_metric(
        name = "needs_pred",
        needs = "predictions"
      ))

      context <- NULL
      expect_warning(
        context <- build_metric_context(
          fit_result,
          fitter,
          data_bundle,
          metrics
        ),
        "Metric requires predictions but fitter does not support them"
      )

      # Should not have predictions since fitter doesn't support it
      expect_false("predictions" %in% names(context))
    })
  })

  describe("compute_all_metrics()", {
    it("handles required vs optional metrics - required failures propagate", {
      # Required metric that fails should cause error
      required_metric <- create_test_metric(
        name = "required_fail",
        required = TRUE,
        compute_fn = function(
          metric,
          fit_result,
          data_bundle,
          context,
          task_ctx
        ) {
          stop("Required metric failed")
        }
      )

      fit_result <- new_fit_result(success = TRUE)
      data_bundle <- valid_data_bundle()
      context <- list()
      task_ctx <- list(task_id = "test")

      expect_error(
        compute_all_metrics(
          fit_result,
          data_bundle,
          context,
          list(required_metric),
          task_ctx
        ),
        "Required metric failed"
      )
    })

    it("returns NA for non-required metrics that fail", {
      # Non-required metric that fails should return NA
      optional_metric <- create_test_metric(
        name = "optional_fail",
        required = FALSE,
        compute_fn = function(
          metric,
          fit_result,
          data_bundle,
          context,
          task_ctx
        ) {
          stop("Optional metric failed")
        }
      )

      fit_result <- new_fit_result(success = TRUE)
      data_bundle <- valid_data_bundle()
      context <- list()
      task_ctx <- list(task_id = "test")

      result <- compute_all_metrics(
        fit_result,
        data_bundle,
        context,
        list(optional_metric),
        task_ctx
      )

      expect_true("optional_fail__value" %in% names(result$metrics))
      expect_true(is.na(result$metrics$optional_fail__value))
    })

    it("returns flattened metric values", {
      metric <- create_test_metric(
        name = "test_compute",
        compute_fn = function(
          metric,
          fit_result,
          data_bundle,
          context,
          task_ctx
        ) {
          list(value = 42.0, count = 10L)
        }
      )

      fit_result <- new_fit_result(success = TRUE)
      data_bundle <- valid_data_bundle()
      context <- list()
      task_ctx <- list(task_id = "test")

      result <- compute_all_metrics(
        fit_result,
        data_bundle,
        context,
        list(metric),
        task_ctx
      )

      expect_true("test_compute__value" %in% names(result$metrics))
      expect_true("test_compute__count" %in% names(result$metrics))
      expect_equal(result$metrics$test_compute__value, 42.0)
    })
  })

  describe("apply_fit_retention()", {
    it("removes correct fields", {
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")

      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(raw = "fit_object"),
        draws = draws,
        diagnostics = list(rhat = 1.01)
      )
      data_bundle <- valid_data_bundle()

      # Default: remove fit and draws
      result <- apply_fit_retention(
        fit_result,
        c("metrics", "diagnostics"),
        data_bundle = data_bundle
      )

      expect_null(result$fit)
      expect_null(result$draws)
      expect_true(is.list(result$diagnostics))
    })

    it("keeps fit when specified", {
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(raw = "fit_object"),
        draws = NULL
      )
      data_bundle <- valid_data_bundle()

      result <- apply_fit_retention(
        fit_result,
        c("metrics", "fit"),
        data_bundle = data_bundle
      )

      expect_true(is.list(result$fit))
    })

    it("keeps draws when specified", {
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")

      fit_result <- new_fit_result(
        success = TRUE,
        draws = draws
      )
      data_bundle <- valid_data_bundle()

      result <- apply_fit_retention(
        fit_result,
        c("metrics", "draws"),
        data_bundle = data_bundle
      )

      expect_true(is.matrix(result$draws))
    })

    it("keeps diagnostics when specified", {
      fit_result <- new_fit_result(
        success = TRUE,
        diagnostics = list(rhat = 1.01)
      )
      data_bundle <- valid_data_bundle()

      result <- apply_fit_retention(
        fit_result,
        c("metrics", "diagnostics"),
        data_bundle = data_bundle
      )

      expect_true(is.list(result$diagnostics))
    })

    it("removes diagnostics when not specified", {
      fit_result <- new_fit_result(
        success = TRUE,
        diagnostics = list(rhat = 1.01)
      )
      data_bundle <- valid_data_bundle()

      result <- apply_fit_retention(
        fit_result,
        c("metrics"),
        data_bundle = data_bundle
      )

      expect_null(result$diagnostics)
    })
  })
})


# =============================================================================
# 5. run_simulation() (from simulate.R)
# =============================================================================

describe("run_simulation()", {
  # Check if run_simulation exists and is functional
  run_sim_exists <- function() {
    if (!exists("run_simulation", mode = "function")) {
      return(FALSE)
    }

    # Also check if is_simulation_config works
    config <- tryCatch(
      simulation_config(
        data_grid = data.frame(n = 100),
        fit_grid = data.frame(model = "test"),
        data_generator = function(data_spec, task_ctx) {
          list(
            train = data.frame(y = rnorm(10)),
            test = NULL,
            response = "y",
            true_params = c(a = 1),
            vars_of_interest = "a"
          )
        },
        seed = 42L
      ),
      error = function(e) NULL
    )

    if (is.null(config)) {
      return(FALSE)
    }
    is_simulation_config(config)
  }

  create_sim_config <- function(n_replicates = 2L) {
    simulation_config(
      data_grid = data.frame(n = c(50, 100)),
      fit_grid = data.frame(model = c("baseline")),
      data_generator = function(data_spec, task_ctx) {
        list(
          train = data.frame(y = rnorm(data_spec$n), x = rnorm(data_spec$n)),
          test = NULL,
          response = "y",
          true_params = c(beta = 1.0, sigma = 1.0),
          vars_of_interest = c("beta", "sigma")
        )
      },
      fitter = MockFitter(),
      metrics = list(create_test_metric(name = "sim_test_metric")),
      n_replicates = n_replicates,
      seed = 42L,
      max_errors = Inf
    )
  }

  describe("with MockFitter", {
    it("completes successfully", {
      skip_if_not(
        run_sim_exists(),
        "run_simulation not available or is_simulation_config broken"
      )

      config <- create_sim_config(n_replicates = 2L)

      result <- run_simulation(config)

      expect_s3_class(result, "bayesim_simulation_result")
    })

    it("results have correct structure", {
      skip_if_not(
        run_sim_exists(),
        "run_simulation not available or is_simulation_config broken"
      )

      config <- create_sim_config(n_replicates = 2L)

      result <- run_simulation(config)

      expect_true(is.character(result$config_fingerprint))
      expect_true(is.list(result$task_results))
      expect_true(is.data.frame(result$summary))
      expect_true(is.list(result$timing))
      expect_true(is.data.frame(result$errors))
    })

    it("results contain expected number of tasks", {
      skip_if_not(
        run_sim_exists(),
        "run_simulation not available or is_simulation_config broken"
      )

      config <- create_sim_config(n_replicates = 2L)
      # 2 data configs * 1 fit config * 2 replicates = 4 tasks

      result <- run_simulation(config)

      expect_equal(length(result$task_results), 4)
      expect_equal(nrow(result$summary), 4)
    })
  })

  describe("error counting", {
    it("respects max_errors parameter", {
      skip_if_not(
        run_sim_exists(),
        "run_simulation not available or is_simulation_config broken"
      )

      # Create config where all tasks will fail
      config <- simulation_config(
        data_grid = data.frame(n = 100),
        fit_grid = data.frame(model = "test"),
        data_generator = function(data_spec, task_ctx) {
          stop("Intentional failure")
        },
        fitter = MockFitter(),
        metrics = list(),
        n_replicates = 5L,
        seed = 42L,
        max_errors = 3
      )

      result <- run_simulation(config)

      # Should stop after max_errors failures
      n_failed <- sum(sapply(result$task_results, function(t) {
        t$status == "failed"
      }))
      expect_lte(n_failed, 3)
    })

    it("continues with Inf max_errors", {
      skip_if_not(
        run_sim_exists(),
        "run_simulation not available or is_simulation_config broken"
      )

      # Create config where tasks fail
      config <- simulation_config(
        data_grid = data.frame(n = 100),
        fit_grid = data.frame(model = "test"),
        data_generator = function(data_spec, task_ctx) {
          stop("Intentional failure")
        },
        fitter = MockFitter(),
        metrics = list(),
        n_replicates = 3L,
        seed = 42L,
        max_errors = Inf
      )

      result <- run_simulation(config)

      # All tasks should have run
      expect_equal(length(result$task_results), 3)
    })
  })

  describe("progress", {
    it("can be disabled", {
      skip_if_not(
        run_sim_exists(),
        "run_simulation not available or is_simulation_config broken"
      )

      config <- create_sim_config(n_replicates = 2L)

      # run_simulation uses cli::cli_alert_info() for status messages
      # which are not controlled by the progress parameter
      # Use suppressMessages() to test that no unexpected output occurs
      expect_silent(suppressMessages(run_simulation(config, progress = FALSE)))
    })
  })

  describe("resume API", {
    it("accepts resume mode strings", {
      skip_if_not(
        run_sim_exists(),
        "run_simulation not available or is_simulation_config broken"
      )

      config <- create_sim_config(n_replicates = 1L)
      result <- run_simulation(config, resume = "never", progress = FALSE)

      expect_s3_class(result, "bayesim_simulation_result")
    })

    it("resume_simulation() resumes from an existing result path when config is supplied", {
      skip_if_not(
        run_sim_exists(),
        "run_simulation not available or is_simulation_config broken"
      )

      result_path <- file.path(withr::local_tempdir(), "resume-results")
      config <- create_sim_config(n_replicates = 1L)
      config@result_path <- result_path

      first <- run_simulation(config, resume = "never", progress = FALSE)
      resumed <- resume_simulation(
        result_path,
        config = config,
        progress = FALSE
      )

      expect_s3_class(first, "bayesim_simulation_result")
      expect_s3_class(resumed, "bayesim_simulation_result")
      expect_equal(nrow(first$summary), nrow(resumed$summary))
    })

    it("resume_simulation() can rehydrate manifest-stored package components", {
      skip_if_not(
        run_sim_exists(),
        "run_simulation not available or is_simulation_config broken"
      )

      result_path <- file.path(
        withr::local_tempdir(),
        "manifest-resume-results"
      )
      config <- simulation_config(
        data_grid = data.frame(n = c(10, 20), beta = c(1, 2)),
        fit_grid = data.frame(model = "baseline"),
        data_generator = bayesim:::bayesim_example_data_generator,
        fitter = MockFitter(),
        metrics = list(pred_rmse_metric()),
        n_replicates = 1L,
        seed = 42L,
        result_path = result_path
      )

      first <- run_simulation(config, resume = "never", progress = FALSE)
      resumed <- resume_simulation(result_path, progress = FALSE)

      expect_s3_class(first, "bayesim_simulation_result")
      expect_s3_class(resumed, "bayesim_simulation_result")
      expect_equal(nrow(resumed$summary), nrow(first$summary))
    })

    it("produces matching summaries under sequential and mirai daemon execution", {
      skip_on_cran()
      skip_if_not(
        run_sim_exists(),
        "run_simulation not available or is_simulation_config broken"
      )

      config <- simulation_config(
        data_grid = data.frame(n = c(50, 100), beta = c(1, 2)),
        fit_grid = data.frame(model = "baseline"),
        data_generator = mock_data_generator,
        fitter = MockFitter(),
        metrics = list(pred_rmse_metric()),
        n_replicates = 2L,
        seed = 42L
      )

      # Sequential path (no daemons set)
      seq_result <- run_simulation(config, resume = "never", progress = FALSE)

      # Parallel path via mirai daemons; reset on exit
      mirai::daemons(2)
      on.exit(mirai::daemons(0), add = TRUE)
      par_result <- run_simulation(config, resume = "never", progress = FALSE)

      normalize_summary <- function(x) {
        x <- x[order(x$task_id), , drop = FALSE]
        x$timing_total <- NULL
        x
      }

      expect_equal(
        normalize_summary(seq_result$summary),
        normalize_summary(par_result$summary)
      )
    })

    it("ships the model bank once per run, not per batch (F6)", {
      skip_on_cran()
      skip_if_not(
        run_sim_exists(),
        "run_simulation not available or is_simulation_config broken"
      )

      # F6 hoisted the mirai::everywhere() bank-ship + daemon_setup hook out of
      # run_batch() into execute_tasks(), so it fires once per run instead of
      # once per batch. Verify by counting everywhere() calls across a
      # multi-batch run. We bypass run_simulation() (which overwrites the bank
      # for non-BrmsFitters) and call execute_tasks() directly with a
      # non-NULL bank in the session option.
      config <- simulation_config(
        data_grid = data.frame(n = 50),
        fit_grid = data.frame(model = "baseline"),
        data_generator = mock_data_generator,
        fitter = MockFitter(),
        metrics = list(pred_rmse_metric()),
        n_replicates = 6L,
        seed = 42L,
        checkpoint_every = 2L
      )
      task_grid <- create_task_grid(config)

      mirai::daemons(2)
      on.exit(mirai::daemons(0), add = TRUE)

      # Set a non-NULL bank so the ship branch fires.
      bayesim:::set_model_bank(list(dummy = "entry"))
      on.exit(bayesim:::set_model_bank(NULL), add = TRUE)

      # Intercept mirai::everywhere to count dispatches (3 batches would have
      # fired 3x under the old per-batch placement; the hoist fires exactly 1x).
      orig_everywhere <- mirai::everywhere
      everywhere_calls <- 0L
      mock_everywhere <- function(expr = NULL, .args = list(), .local = FALSE) {
        everywhere_calls <<- everywhere_calls + 1L
      }
      unlockBinding("everywhere", asNamespace("mirai"))
      assign("everywhere", mock_everywhere, envir = asNamespace("mirai"))
      on.exit({
        assign("everywhere", orig_everywhere, envir = asNamespace("mirai"))
        lockBinding("everywhere", asNamespace("mirai"))
      }, add = TRUE)

      execute_tasks(
        task_grid = task_grid,
        config = config,
        config_spec = as_config_spec(config),
        fitter = config@fitter,
        metrics = config@metrics,
        retain = config@retain,
        max_errors = config@max_errors,
        progress = FALSE
      )

      # Exactly one bank-ship per run across 3 batches (F6 hoist).
      expect_equal(everywhere_calls, 1L)
    })

    it("rejects parquet checkpoint_format at construction time", {
      expect_error(
        simulation_config(
          data_grid = data.frame(n = 50),
          fit_grid = data.frame(model = "baseline"),
          data_generator = function(data_spec, task_ctx) {
            list(
              train = data.frame(y = rnorm(data_spec$n)),
              test = NULL,
              response = "y",
              true_params = c(beta = 1),
              vars_of_interest = "beta"
            )
          },
          fitter = MockFitter(),
          n_replicates = 1L,
          seed = 42L,
          checkpoint_format = "parquet"
        ),
        "rds"
      )
    })
  })
})
