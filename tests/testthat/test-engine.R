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
# 1. RNG Management (from rng.R)
# =============================================================================

describe("RNG Management", {
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

    it("restores every RNG kind component and an absent seed", {
      original_kind <- RNGkind()
      original_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
        get(".Random.seed", envir = .GlobalEnv)
      } else {
        NULL
      }
      on.exit(
        {
          do.call(base::RNGkind, as.list(original_kind))
          if (is.null(original_seed)) {
            if (exists(".Random.seed", envir = .GlobalEnv)) {
              rm(".Random.seed", envir = .GlobalEnv)
            }
          } else {
            assign(".Random.seed", original_seed, envir = .GlobalEnv)
          }
        },
        add = TRUE
      )

      RNGkind("Mersenne-Twister", "Inversion", "Rejection")
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
      before_kind <- RNGkind()

      create_task_rng_streams(42, 3)

      expect_identical(RNGkind(), before_kind)
      expect_false(exists(".Random.seed", envir = .GlobalEnv))
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

  describe("get_task_count_summary()", {
    it("returns aligned zero-filled counts for every status", {
      grid <- data.frame(status = c("success", "success", "failed"))

      counts <- get_task_count_summary(grid)

      expect_equal(
        counts,
        c(pending = 0L, success = 2L, failed = 1L, skipped = 0L)
      )
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

      task_grid$status[1:3] <- c("success", "failed", "success")
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

      task_grid$status[1:3] <- c("success", "failed", "success")
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
      expect_s7_object(posterior_summary_metric())
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

    it("stays total when error-capture helpers themselves fail", {
      # Mirai daemons without an installed bayesim (source-loaded controller)
      # can fail to resolve the handler's package helpers; the handler must
      # fall back to base-R capture instead of dying (PR #49 review).
      task <- create_test_task()
      config_spec <- list(
        data_generator = function(data_spec, task_ctx) {
          stop("Data generation error")
        }
      )

      testthat::local_mocked_bindings(
        capture_error_info = function(e) stop("namespace unavailable"),
        .package = "bayesim"
      )
      result <- run_task_safe(task, config_spec, MockFitter(), list())

      expect_s3_class(result, "bayesim_task_result")
      expect_equal(result$status, "failed")
      expect_equal(result$task_id, "d001_f001_r00001")
      # run_task() itself calls capture_error_info() on generator errors, so
      # the captured condition may be the mock's error rather than the
      # generator's; either way the fallback must record a non-empty message
      # and flag the degradation.
      expect_type(result$error$error_message, "character")
      expect_true(nzchar(result$error$error_message))
      expect_match(result$error$handler_error, "fell back to base capture")
      expect_false(result$error$fatal)
      expect_type(result$error$condition_class, "character")
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

    it("warns once when a metric needs epred the fitter does not support", {
      # MockFitter supports loo but inherits supports_epred = FALSE — the
      # silent-NA shape from #62: r2_loo/rmse_loo must explain their
      # degradation instead of returning all-NA columns.
      bayesim:::.reset_warn_once()
      fitter <- MockFitter()
      data_bundle <- valid_data_bundle()

      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )

      context <- NULL
      expect_warning(
        context <- build_metric_context(
          fit_result,
          fitter,
          data_bundle,
          list(r2_loo_metric())
        ),
        "Metric requires epred but fitter does not support it"
      )
      # The loo summary is still built (loo is supported); only the epred part
      # of the LOO context is absent.
      expect_false(is.null(context$loo))
      expect_true(is.null(context$loo_epred))

      # Warn once per run, not per task.
      expect_silent(
        build_metric_context(
          fit_result,
          fitter,
          data_bundle,
          list(r2_loo_metric())
        )
      )
    })

    it("warns once when an epred-supporting fitter fails to produce epred", {
      # supports_epred = TRUE but predict_epred() errors at run time: the
      # unsupported-capability warning must not fire; the LOO-context seam
      # explains the NA degradation instead (#62). S7 keeps the parent's
      # property default when a subclass redeclares it, so set the flag on
      # the instance.
      bayesim:::.reset_warn_once()
      BrokenEpredFitter <- S7::new_class(
        "BrokenEpredFitter",
        parent = MockFitter
      )
      S7::method(predict_epred, BrokenEpredFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        stop("boom")
      }
      fitter <- BrokenEpredFitter()
      fitter@supports_epred <- TRUE

      data_bundle <- valid_data_bundle()
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )

      context <- NULL
      expect_warning(
        context <- build_metric_context(
          fit_result,
          fitter,
          data_bundle,
          list(r2_loo_metric())
        ),
        "epred matrix unavailable"
      )
      expect_false(is.null(context$loo))
      expect_true(is.null(context$loo_epred))
    })

    it("degrades a wrong-shaped predict_epred return to the warn-once path", {
      # A vector or dimension-mismatched matrix would otherwise bypass the
      # seam warning and die as a generic metric error inside loo::E_loo().
      bayesim:::.reset_warn_once()
      BadShapeEpredFitter <- S7::new_class(
        "BadShapeEpredFitter",
        parent = MockFitter
      )
      S7::method(predict_epred, BadShapeEpredFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        matrix(rnorm(10), nrow = 5, ncol = 2)
      }
      fitter <- BadShapeEpredFitter()
      fitter@supports_epred <- TRUE

      data_bundle <- valid_data_bundle()
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )

      context <- NULL
      expect_warning(
        context <- build_metric_context(
          fit_result,
          fitter,
          data_bundle,
          list(r2_loo_metric())
        ),
        "epred matrix unavailable"
      )
      # The malformed matrix must not reach the metric context.
      expect_true(is.null(context$loo_epred))
      expect_false(is.null(context$loo))
    })

    it("builds epred for a metric declaring needs = \"epred\" alone (#68)", {
      # epred used to be built only inside the LOO branch, so an epred-only
      # metric silently received context$loo_epred = NULL. The matrix must
      # now be built directly, without the LOO machinery.
      bayesim:::.reset_warn_once()
      EpredCapableFitter <- S7::new_class(
        "EpredCapableFitter",
        parent = MockFitter
      )
      S7::method(predict_epred, EpredCapableFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        matrix(rnorm(500), nrow = 50, ncol = 10)
      }
      fitter <- EpredCapableFitter()
      fitter@supports_epred <- TRUE

      data_bundle <- valid_data_bundle()
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )
      # Build the metric outside expect_silent: registering a fresh S7 class
      # emits messages that expect_silent would flag.
      metric <- list(create_test_metric(name = "epred_only", needs = "epred"))

      context <- expect_silent(build_metric_context(
        fit_result,
        fitter,
        data_bundle,
        metric
      ))

      # epred on the training set: 50 draws x 10 observations.
      expect_equal(dim(context$loo_epred), c(50, 10))
      # No "loo" need: the LOO context must not be built. Exact [[ ]] access —
      # $loo would partially match loo_epred.
      expect_true(is.null(context[["loo"]]))
      expect_true(is.null(context[["loo_psis"]]))
      expect_true(is.null(context[["loo_psis_ll"]]))
    })

    it("skips the PSIS machinery when only \"loo\" is needed (#69)", {
      # needs = "loo" alone (the elpd_loo configuration) must pay for the
      # loo_fit() summary only: the train-set log-lik matrix, r_eff, the PSIS
      # object, and epred exist solely to feed the weighted-prediction
      # machinery, which no such metric reads.
      CountingLooFitter <- S7::new_class(
        "CountingLooFitter",
        parent = MockFitter
      )
      calls <- new.env(parent = emptyenv())
      calls$log_lik <- 0L
      calls$epred <- 0L
      S7::method(log_lik_matrix, CountingLooFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        calls$log_lik <- calls$log_lik + 1L
        matrix(rnorm(500), nrow = 50, ncol = 10)
      }
      S7::method(predict_epred, CountingLooFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        calls$epred <- calls$epred + 1L
        matrix(rnorm(500), nrow = 50, ncol = 10)
      }
      fitter <- CountingLooFitter()
      fitter@supports_epred <- TRUE

      data_bundle <- valid_data_bundle()
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )

      context <- expect_silent(build_metric_context(
        fit_result,
        fitter,
        data_bundle,
        list(elpd_loo_metric())
      ))

      # The summary is delivered (MockFitter: elpd = -n_obs * 2)...
      expect_false(is.null(context[["loo"]]))
      expect_equal(context$loo$elpd, -20)
      # ...the PSIS machinery is not computed, not merely absent from the
      # context.
      expect_null(context$loo_psis)
      expect_null(context$loo_psis_ll)
      expect_null(context$loo_epred)
      expect_equal(calls$log_lik, 0L)
      expect_equal(calls$epred, 0L)
    })

    it("builds the full PSIS machinery when \"epred\" is declared alongside \"loo\" (#69)", {
      # rmse_loo/r2_loo declare needs = c("loo", "epred"): the weighted-
      # prediction machinery must still be computed and shared.
      CountingLooFitter2 <- S7::new_class(
        "CountingLooFitter2",
        parent = MockFitter
      )
      calls2 <- new.env(parent = emptyenv())
      calls2$log_lik <- 0L
      calls2$epred <- 0L
      S7::method(log_lik_matrix, CountingLooFitter2) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        calls2$log_lik <- calls2$log_lik + 1L
        matrix(rnorm(500), nrow = 50, ncol = 10)
      }
      S7::method(predict_epred, CountingLooFitter2) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        calls2$epred <- calls2$epred + 1L
        matrix(rnorm(500), nrow = 50, ncol = 10)
      }
      fitter <- CountingLooFitter2()
      fitter@supports_epred <- TRUE

      data_bundle <- valid_data_bundle()
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )

      context <- expect_silent(build_metric_context(
        fit_result,
        fitter,
        data_bundle,
        list(rmse_loo_metric())
      ))

      expect_false(is.null(context[["loo"]]))
      expect_true(inherits(context$loo_psis, "psis"))
      expect_true(is.matrix(context$loo_psis_ll))
      expect_equal(dim(context$loo_epred), c(50, 10))
      expect_equal(calls2$log_lik, 1L)
      expect_equal(calls2$epred, 1L)
    })

    it("computes the train-set log-lik matrix once on the PSIS path (#73)", {
      # loo_fit() used to recompute the train-set log-lik matrix (and r_eff)
      # that build_loo_context() had just computed. The matrix now travels to
      # loo_fit() via its log_lik argument. The counting fitter mirrors the
      # CmdStanFitter/LinearRegressionFitter shape: their loo_fit() derives
      # from log_lik_matrix() when called without a precomputed matrix
      # (BrmsFitter instead falls back to brms' stored model frame).
      SharingLooFitter <- S7::new_class("SharingLooFitter", parent = MockFitter)
      calls <- new.env(parent = emptyenv())
      calls$log_lik <- 0L
      S7::method(log_lik_matrix, SharingLooFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        calls$log_lik <- calls$log_lik + 1L
        # Deterministic matrix so the handoff can be asserted by identity.
        calls$last_ll <- matrix(seq_len(500) / 500, nrow = 50, ncol = 10)
      }
      S7::method(loo_fit, SharingLooFitter) <- function(
        fitter,
        fit_result,
        log_lik = NULL,
        save_psis = FALSE
      ) {
        ll <- log_lik %||% log_lik_matrix(fitter, fit_result)
        calls$received_passed_ll <- !is.null(log_lik)
        calls$received_save_psis <- save_psis
        list(
          elpd = -20,
          p_loo = 3,
          elpd_se = 1.5,
          pareto_k = runif(ncol(ll)),
          r_eff = NULL
        )
      }
      S7::method(predict_epred, SharingLooFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        matrix(rnorm(500), nrow = 50, ncol = 10)
      }
      fitter <- SharingLooFitter()
      fitter@supports_epred <- TRUE

      data_bundle <- valid_data_bundle()
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )

      context <- expect_silent(build_metric_context(
        fit_result,
        fitter,
        data_bundle,
        list(rmse_loo_metric())
      ))

      # Exactly one train-set log-lik computation feeds both the loo_fit()
      # summary and the PSIS object...
      expect_equal(calls$log_lik, 1L)
      # ...and loo_fit() received the matrix the context computed.
      expect_true(isTRUE(calls$received_passed_ll))
      # The PSIS path asks loo_fit() to retain its PSIS object (#76); this
      # fitter returns none, so the loo::psis() fallback built the context's
      # object.
      expect_true(isTRUE(calls$received_save_psis))
      expect_identical(context$loo_psis_ll, calls$last_ll)
      expect_false(is.null(context[["loo"]]))
      expect_true(inherits(context$loo_psis, "psis"))
    })

    it("reuses loo_fit()'s psis_object instead of re-smoothing the tails (#76)", {
      # The PSIS tail smoothing used to run twice: once inside loo_fit()'s
      # loo::loo() and once more in build_loo_context()'s loo::psis(). Now
      # loo_fit() receives save_psis = TRUE on the PSIS path and its returned
      # psis_object — identical to loo::psis(-ll, r_eff) on the same matrix —
      # is reused directly. This fitter mirrors the built-in implementations:
      # it runs the real loo machinery and returns the retained object.
      PsisSavingLooFitter <- S7::new_class(
        "PsisSavingLooFitter",
        parent = MockFitter
      )
      calls <- new.env(parent = emptyenv())
      S7::method(log_lik_matrix, PsisSavingLooFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        # Deterministic matrix so the summary and the reused PSIS object are
        # fitted from a known input.
        matrix(seq_len(500) / 500, nrow = 50, ncol = 10)
      }
      S7::method(loo_fit, PsisSavingLooFitter) <- function(
        fitter,
        fit_result,
        log_lik = NULL,
        save_psis = FALSE
      ) {
        ll <- log_lik %||% log_lik_matrix(fitter, fit_result)
        calls$received_save_psis <- save_psis
        loo_result <- suppressWarnings(loo::loo(ll, save_psis = save_psis))
        calls$returned_psis <- loo_result$psis_object
        calls$elpd <- loo_result$estimates["elpd_loo", "Estimate"]
        list(
          elpd = loo_result$estimates["elpd_loo", "Estimate"],
          p_loo = loo_result$estimates["p_loo", "Estimate"],
          elpd_se = loo_result$estimates["elpd_loo", "SE"],
          pareto_k = loo::pareto_k_values(loo_result),
          r_eff = NULL,
          psis_object = loo_result$psis_object
        )
      }
      S7::method(predict_epred, PsisSavingLooFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        matrix(rnorm(500), nrow = 50, ncol = 10)
      }
      fitter <- PsisSavingLooFitter()
      fitter@supports_epred <- TRUE

      data_bundle <- valid_data_bundle()
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )

      context <- expect_silent(build_metric_context(
        fit_result,
        fitter,
        data_bundle,
        list(rmse_loo_metric())
      ))

      # The engine asked loo_fit() to retain its PSIS object...
      expect_true(isTRUE(calls$received_save_psis))
      # ...and reused exactly that object (identity, not merely equality: a
      # re-smoothing loo::psis() run would allocate a distinct object).
      expect_identical(context$loo_psis, calls$returned_psis)
      expect_true(inherits(context$loo_psis, "psis"))
      expect_false(is.null(context[["loo"]]))
      expect_equal(context$loo$elpd, calls$elpd)

      # Standalone calls keep the default and get no PSIS object back.
      standalone <- loo_fit(fitter, fit_result)
      expect_null(standalone$psis_object)
      expect_false(isTRUE(calls$received_save_psis))
    })

    it("falls back to loo::psis() when loo_fit() returns no psis_object (#76)", {
      # SharingLooFitter (the #73 test above) returns no psis_object, like a
      # custom fitter written before #76; build_metric_context still built a
      # usable PSIS object there via the direct loo::psis() route. Pin the
      # same fallback for a fitter whose psis_object does not match the
      # context's log-lik matrix (e.g. it ignored the supplied log_lik).
      MismatchedPsisFitter <- S7::new_class(
        "MismatchedPsisFitter",
        parent = MockFitter
      )
      S7::method(log_lik_matrix, MismatchedPsisFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        matrix(seq_len(500) / 500, nrow = 50, ncol = 10)
      }
      S7::method(loo_fit, MismatchedPsisFitter) <- function(
        fitter,
        fit_result,
        log_lik = NULL,
        save_psis = FALSE
      ) {
        ll <- log_lik %||% log_lik_matrix(fitter, fit_result)
        # PSIS object fitted from a same-width but differently-sized matrix
        # (40 draws vs the context's 50, same 10 observations): only the
        # draws dimension betrays it, and the engine must reject it and
        # re-fit from its own matrix instead of erroring opaquely inside
        # loo::E_loo().
        other <- matrix(rnorm(400, -1, 0.5), nrow = 40, ncol = 10)
        bad_psis <- suppressWarnings(loo::loo(
          other,
          save_psis = TRUE
        ))$psis_object
        list(
          elpd = -20,
          p_loo = 3,
          elpd_se = 1.5,
          pareto_k = runif(ncol(ll)),
          r_eff = NULL,
          psis_object = bad_psis
        )
      }
      S7::method(predict_epred, MismatchedPsisFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        matrix(rnorm(500), nrow = 50, ncol = 10)
      }
      fitter <- MismatchedPsisFitter()
      fitter@supports_epred <- TRUE

      data_bundle <- valid_data_bundle()
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )

      context <- expect_silent(build_metric_context(
        fit_result,
        fitter,
        data_bundle,
        list(rmse_loo_metric())
      ))

      # A usable PSIS object fitted from the context's own 50 x 10 matrix
      # (the mismatched 40-draw object was rejected).
      expect_true(inherits(context$loo_psis, "psis"))
      expect_identical(dim(context$loo_psis$log_weights), c(50L, 10L))
    })

    it("builds epred without LOO support when \"loo\" is also declared (#68)", {
      # needs = c("loo", "epred") on a supports_loo = FALSE fitter: the LOO
      # branch is skipped (warn-once), but epred no longer disappears with it.
      bayesim:::.reset_warn_once()
      EpredNoLooFitter <- S7::new_class(
        "EpredNoLooFitter",
        parent = MockFitter
      )
      S7::method(predict_epred, EpredNoLooFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        matrix(rnorm(500), nrow = 50, ncol = 10)
      }
      fitter <- EpredNoLooFitter()
      fitter@supports_epred <- TRUE
      fitter@supports_loo <- FALSE

      data_bundle <- valid_data_bundle()
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )

      context <- NULL
      expect_warning(
        context <- build_metric_context(
          fit_result,
          fitter,
          data_bundle,
          list(create_test_metric(name = "both", needs = c("loo", "epred")))
        ),
        "Metric requires loo but fitter does not support it"
      )
      expect_true(is.null(context[["loo"]]))
      expect_equal(dim(context$loo_epred), c(50, 10))
    })

    it("elpd_loo composed with an epred metric NA-degrades without LOO support (#68)", {
      # Regression for the $-partial-match hazard (pullfrog review on #70):
      # this branch can produce a context carrying loo_epred without a loo
      # summary, and a $-reading metric must not receive the epred matrix
      # via context$loo.
      bayesim:::.reset_warn_once()
      EpredNoLooComposeFitter <- S7::new_class(
        "EpredNoLooComposeFitter",
        parent = MockFitter
      )
      S7::method(predict_epred, EpredNoLooComposeFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        matrix(rnorm(500), nrow = 50, ncol = 10)
      }
      fitter <- EpredNoLooComposeFitter()
      fitter@supports_epred <- TRUE
      fitter@supports_loo <- FALSE

      data_bundle <- valid_data_bundle()
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )

      metrics <- list(elpd_loo_metric(), rmse_loo_metric())
      context <- NULL
      expect_warning(
        context <- build_metric_context(
          fit_result,
          fitter,
          data_bundle,
          metrics
        ),
        "Metric requires loo but fitter does not support it"
      )
      # The pinned NULL binding must win over $ partial matching of loo_epred.
      expect_null(context$loo)
      expect_equal(dim(context$loo_epred), c(50, 10))

      # The loo-only metric degrades to clean NAs, not an error from being
      # handed the epred matrix.
      out <- compute_metric(
        elpd_loo_metric(),
        fit_result,
        data_bundle,
        context,
        list(task_id = "t1", rep_idx = 1L)
      )
      expect_true(all(is.na(out)))
    })

    it("builds epred when the LOO context bails at the train-set log-lik (#68)", {
      # build_loo_context() returns a partial context (loo summary only) when
      # log_lik_matrix() fails, never attempting epred there; the direct
      # fallback must still deliver the matrix (coderabbit review on #70).
      bayesim:::.reset_warn_once()
      BrokenLlFitter <- S7::new_class("BrokenLlFitter", parent = MockFitter)
      S7::method(log_lik_matrix, BrokenLlFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        stop("ll boom")
      }
      S7::method(predict_epred, BrokenLlFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        matrix(rnorm(500), nrow = 50, ncol = 10)
      }
      fitter <- BrokenLlFitter()
      fitter@supports_epred <- TRUE

      data_bundle <- valid_data_bundle()
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )
      metric <- list(create_test_metric(
        name = "both",
        needs = c("loo", "epred")
      ))

      context <- expect_silent(build_metric_context(
        fit_result,
        fitter,
        data_bundle,
        metric
      ))

      # Partial LOO context: the summary is present, the PSIS machinery is
      # not (rmse_loo/r2_loo NA-degrade on the missing psis/log_lik)...
      expect_false(is.null(context[["loo"]]))
      expect_null(context$loo_psis)
      expect_null(context$loo_psis_ll)
      # ...but epred does not depend on the log-lik and must be delivered.
      expect_equal(dim(context$loo_epred), c(50, 10))
    })

    it("warns once when an epred-only metric's predict_epred fails (#68)", {
      bayesim:::.reset_warn_once()
      BrokenEpredOnlyFitter <- S7::new_class(
        "BrokenEpredOnlyFitter",
        parent = MockFitter
      )
      S7::method(predict_epred, BrokenEpredOnlyFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        stop("boom")
      }
      fitter <- BrokenEpredOnlyFitter()
      fitter@supports_epred <- TRUE

      data_bundle <- valid_data_bundle()
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )

      metric <- list(create_test_metric(name = "epred_only", needs = "epred"))

      context <- NULL
      expect_warning(
        context <- build_metric_context(
          fit_result,
          fitter,
          data_bundle,
          metric
        ),
        "predict_epred\\(\\) failed"
      )
      expect_true(is.null(context$loo_epred))
      # Warn once per run, not per task.
      expect_silent(build_metric_context(
        fit_result,
        fitter,
        data_bundle,
        metric
      ))
    })

    it("degrades a wrong-shaped epred through the warn-once path without loo (#68)", {
      # Outside the LOO context there is no log-lik matrix to align draws
      # against, but the observation count must still match the training set.
      bayesim:::.reset_warn_once()
      BadShapeEpredOnlyFitter <- S7::new_class(
        "BadShapeEpredOnlyFitter",
        parent = MockFitter
      )
      S7::method(predict_epred, BadShapeEpredOnlyFitter) <- function(
        fitter,
        fit_result,
        newdata = NULL
      ) {
        matrix(rnorm(10), nrow = 5, ncol = 2)
      }
      fitter <- BadShapeEpredOnlyFitter()
      fitter@supports_epred <- TRUE

      data_bundle <- valid_data_bundle()
      draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
      colnames(draws) <- c("alpha", "beta")
      fit_result <- new_fit_result(
        success = TRUE,
        fit = list(data_bundle = data_bundle, seed = 42L, n_obs = 10),
        draws = draws
      )

      context <- NULL
      expect_warning(
        context <- build_metric_context(
          fit_result,
          fitter,
          data_bundle,
          list(create_test_metric(name = "epred_only", needs = "epred"))
        ),
        "predict_epred\\(\\) failed"
      )
      # The malformed matrix must not reach the metric context.
      expect_true(is.null(context$loo_epred))
      expect_true(is.null(context[["loo"]]))
    })
  })

  describe("compute_all_metrics()", {
    it("isolates stochastic metrics from metric ordering", {
      StochasticTestMetric <- S7::new_class(
        "StochasticTestMetric",
        parent = Metric,
        properties = list(
          name = S7::new_property(S7::class_character),
          needs = S7::new_property(S7::class_character, default = character()),
          required = S7::new_property(S7::class_logical, default = FALSE)
        )
      )
      S7::method(compute_metric, StochasticTestMetric) <- function(
        metric,
        fit_result,
        data_bundle,
        context,
        task_ctx
      ) {
        list(value = stats::runif(1))
      }
      metric_a <- StochasticTestMetric(name = "random_a")
      metric_b <- StochasticTestMetric(name = "random_b")
      args <- list(
        fit_result = new_fit_result(success = TRUE),
        data_bundle = valid_data_bundle(),
        context = list(),
        task_ctx = list(task_id = "t1", seed = 123L)
      )

      forward <- do.call(
        compute_all_metrics,
        c(args, list(metrics = list(metric_a, metric_b)))
      )
      reverse <- do.call(
        compute_all_metrics,
        c(args, list(metrics = list(metric_b, metric_a)))
      )

      expect_identical(
        forward$metrics[c("random_a__value", "random_b__value")],
        reverse$metrics[c("random_a__value", "random_b__value")]
      )
    })

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

  describe("prediction retention", {
    it("copies metric-context predictions onto the task result", {
      task_result <- new_task_result(
        task_id = "t1",
        status = "success",
        metrics = list(value = 1),
        timing = list(total = 0)
      )
      fit_result <- new_fit_result(success = TRUE)
      predictions <- list(
        predicted_mean = c(1, 2),
        predicted_samples = matrix(1:4, nrow = 2)
      )

      retained <- apply_task_retention(
        task_result,
        fit_result,
        valid_data_bundle(),
        retain = c("metrics", "predictions"),
        predictions = predictions
      )

      expect_identical(retained$predictions, predictions)
    })

    it("computes retained predictions even when no metric requests them", {
      config <- simulation_config(
        data_grid = data.frame(n = 10L),
        fit_grid = data.frame(model = "linear"),
        data_generator = function(data_spec, task_ctx) {
          x <- stats::rnorm(data_spec$n)
          list(
            train = data.frame(y = 1 + x, x = x),
            test = NULL,
            response = "y",
            true_params = c(Intercept = 1, x = 1, sigma = 0.1),
            vars_of_interest = c("Intercept", "x", "sigma")
          )
        },
        fitter = LinearRegressionFitter(n_draws = 10L),
        metrics = list(),
        n_replicates = 1L,
        seed = 42L,
        retain = c("metrics", "predictions")
      )

      result <- run_simulation(config, resume = "never", progress = FALSE)

      expect_true(is.list(result$task_results[[1]]$predictions))
      expect_length(result$task_results[[1]]$predictions$predicted_mean, 10L)
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
      on.exit(
        {
          assign("everywhere", orig_everywhere, envir = asNamespace("mirai"))
          lockBinding("everywhere", asNamespace("mirai"))
        },
        add = TRUE
      )

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
