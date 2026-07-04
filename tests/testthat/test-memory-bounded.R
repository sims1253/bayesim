# test-memory-bounded.R
# Tests for memory-bounded execution functionality
# testthat 3e style with describe/it blocks

library(bayesim)

# =============================================================================
# Helper functions
# =============================================================================

#' Create a simple test task result with optional heavy objects
create_heavy_task_result <- function(
  task_id = "t1",
  status = "success",
  include_fit = FALSE,
  include_draws = FALSE,
  include_data = FALSE
) {
  result <- list(
    task_id = task_id,
    status = status,
    metrics = list(rmse = 0.05, bias = 0.01),
    diagnostics = list(rhat = 1.01, ess = 1000),
    timing = list(total = 1.5),
    error = NULL,
    warnings = c("test warning")
  )

  if (include_fit) {
    result$fit <- list(heavy_object = rnorm(10000))
  }
  if (include_draws) {
    result$draws <- matrix(rnorm(4000), nrow = 1000, ncol = 4)
  }
  if (include_data) {
    result$data <- list(train = data.frame(x = 1:100, y = rnorm(100)))
  }

  result
}

# =============================================================================
# Tests for lighten_task_result()
# =============================================================================

describe("lighten_task_result()", {
  it("preserves essential fields (task_id, status, metrics, timing, error)", {
    heavy <- create_heavy_task_result(
      task_id = "test_task",
      status = "success",
      include_fit = TRUE,
      include_draws = TRUE
    )

    light <- lighten_task_result(heavy, retain = c("metrics"))

    expect_equal(light$task_id, "test_task")
    expect_equal(light$status, "success")
    expect_equal(light$metrics$rmse, 0.05)
    expect_equal(light$timing$total, 1.5)
    expect_null(light$error)
  })

  it("removes heavy objects (fit, draws, data) regardless of retain", {
    heavy <- create_heavy_task_result(
      include_fit = TRUE,
      include_draws = TRUE,
      include_data = TRUE
    )

    # Even with "debug" retention, lighten removes heavy objects
    light <- lighten_task_result(
      heavy,
      retain = c("metrics", "diagnostics", "fit", "draws", "data")
    )

    expect_null(light$fit)
    expect_null(light$draws)
    expect_null(light$data)
  })

  it("keeps diagnostics when 'diagnostics' is in retain", {
    heavy <- create_heavy_task_result()

    light <- lighten_task_result(heavy, retain = c("metrics", "diagnostics"))

    expect_equal(light$diagnostics$rhat, 1.01)
    expect_equal(light$diagnostics$ess, 1000)
  })

  it("removes diagnostics when 'diagnostics' is NOT in retain", {
    heavy <- create_heavy_task_result()

    light <- lighten_task_result(heavy, retain = c("metrics"))

    expect_null(light$diagnostics)
  })

  it("keeps warnings when 'warnings' is in retain", {
    heavy <- create_heavy_task_result()

    light <- lighten_task_result(heavy, retain = c("metrics", "warnings"))

    expect_equal(light$warnings, "test warning")
  })

  it("removes warnings when 'warnings' is NOT in retain", {
    heavy <- create_heavy_task_result()

    light <- lighten_task_result(heavy, retain = c("metrics"))

    expect_null(light$warnings)
  })

  it("returns NULL for NULL input", {
    expect_null(lighten_task_result(NULL, retain = c("metrics")))
  })

  it("produces smaller object than input", {
    heavy <- create_heavy_task_result(
      include_fit = TRUE,
      include_draws = TRUE
    )

    light <- lighten_task_result(heavy, retain = c("metrics", "diagnostics"))

    heavy_size <- estimate_size(heavy)
    light_size <- estimate_size(light)

    expect_true(light_size < heavy_size)
  })
})

# =============================================================================
# Tests for retention profile updates
# =============================================================================

describe("Retention profiles", {
  it("standard profile now includes warnings", {
    standard <- resolve_retention("standard")

    expect_true("metrics" %in% standard)
    expect_true("diagnostics" %in% standard)
    expect_true("warnings" %in% standard)
  })

  it("minimal profile only includes metrics", {
    minimal <- resolve_retention("minimal")

    expect_equal(minimal, c("metrics"))
  })

  it("debug profile includes all options", {
    debug <- resolve_retention("debug")

    expect_true("metrics" %in% debug)
    expect_true("diagnostics" %in% debug)
    expect_true("draws" %in% debug)
    expect_true("predictions" %in% debug)
    expect_true("fit" %in% debug)
    expect_true("data" %in% debug)
    expect_true("warnings" %in% debug)
  })
})

# =============================================================================
# Tests for checkpoint_every parameter (B4: chunk_size merged into it)
# =============================================================================

describe("simulation_config() checkpoint_every parameter", {
  .gen <- function(data_spec, task_ctx) {
    list(
      train = data.frame(x = 1:10, y = 1:10),
      test = NULL,
      response = "y",
      true_params = c(a = 1),
      vars_of_interest = "a",
      meta = list()
    )
  }

  it("accepts checkpoint_every parameter", {
    config <- simulation_config(
      data_grid = data.frame(n = 100),
      fit_grid = data.frame(model = "test"),
      data_generator = .gen,
      fitter = NULL,
      n_replicates = 2L,
      seed = 42L,
      checkpoint_every = 10L
    )
    expect_equal(config@checkpoint_every, 10L)
  })

  it("B4: chunk_size and max_in_memory are no longer arguments", {
    # Pre-release API: the merged knobs were removed entirely (no shim).
    expect_error(
      simulation_config(
        data_grid = data.frame(n = 100),
        fit_grid = data.frame(model = "test"),
        data_generator = .gen,
        fitter = NULL,
        n_replicates = 2L,
        seed = 42L,
        chunk_size = 10L
      ),
      class = "simpleError"
    )
  })

  it("validates checkpoint_every is positive integer", {
    expect_error(
      simulation_config(
        data_grid = data.frame(n = 100),
        fit_grid = data.frame(model = "test"),
        data_generator = .gen,
        fitter = NULL,
        n_replicates = 2L,
        seed = 42L,
        checkpoint_every = 0L
      ),
      "checkpoint_every must be a positive integer"
    )
  })
})
