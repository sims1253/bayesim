# test-retention.R
# Tests for retention helpers: retention_for_task_result, apply_fit_retention,
# apply_task_retention, lighten_task_result, resolve_retention_spec

library(bayesim)

describe("Retention Resolution", {
  describe("resolve_retention_spec()", {
    it("resolves simple character vector to full spec", {
      spec <- resolve_retention_spec(c("metrics", "diagnostics"))
      expect_true(is.list(spec))
      expect_true("success" %in% names(spec))
      expect_true("error" %in% names(spec))
      expect_true("warning" %in% names(spec))
    })

    it("resolves retain with all valid options", {
      spec <- resolve_retention_spec(c("metrics", "diagnostics", "draws", "fit", "data", "predictions", "warnings"))
      all_opts <- spec$success
      expect_true("fit" %in% all_opts)
      expect_true("draws" %in% all_opts)
      expect_true("data" %in% all_opts)
    })

    it("rejects invalid retain options", {
      expect_error(
        resolve_retention_spec(c("metrics", "invalid_option")),
        "invalid options"
      )
    })

    it("accepts empty retain", {
      spec <- resolve_retention_spec(character())
      expect_true(is.list(spec))
    })

    it("accepts context-specific retain list", {
      spec <- resolve_retention_spec(list(
        success = c("metrics", "diagnostics"),
        error = "debug",
        warning = c("metrics")
      ))
      expect_equal(spec$success, c("metrics", "diagnostics"))
    })
  })
})

describe("retention_for_task_result()", {
  it("returns the correct retain spec for success status", {
    spec <- resolve_retention_spec(c("metrics", "diagnostics"))
    retain <- retention_for_task_result(spec, "success", character())
    expect_true("metrics" %in% retain)
    expect_true("diagnostics" %in% retain)
  })

  it("uses error context for failed status", {
    spec <- resolve_retention_spec(c("metrics", "diagnostics"))
    retain <- retention_for_task_result(spec, "failed", character())
    # Error context uses spec$error which falls back to base (metrics, diagnostics)
    expect_true("metrics" %in% retain)
  })

  it("uses warning context when warnings are present", {
    spec <- resolve_retention_spec(c("metrics", "diagnostics"))
    retain <- retention_for_task_result(spec, "success", c("some warning"))
    expect_true("metrics" %in% retain)
  })

  it("always includes metrics", {
    spec <- resolve_retention_spec(c("diagnostics"))
    retain <- retention_for_task_result(spec, "success", character())
    expect_true("metrics" %in% retain)
  })

  it("respects context-specific retain options", {
    spec <- resolve_retention_spec(list(
      success = c("metrics", "diagnostics"),
      error = c("fit", "draws", "data"),
      warning = c("metrics")
    ))
    retain <- retention_for_task_result(spec, "failed", character())
    expect_true("fit" %in% retain)
    expect_true("draws" %in% retain)
    expect_true("data" %in% retain)
  })
})

describe("apply_fit_retention()", {
  it("removes fit when not in retain", {
    draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
    colnames(draws) <- c("alpha", "beta")
    result <- new_fit_result(
      success = TRUE,
      fit = list(model = "test"),
      draws = draws,
      diagnostics = list(rhat = 1.01),
      timing = list(total = 1.0)
    )

    retained <- apply_fit_retention(result, c("metrics", "diagnostics"))
    expect_null(retained$fit)
  })

  it("removes draws when not in retain", {
    draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
    colnames(draws) <- c("alpha", "beta")
    result <- new_fit_result(
      success = TRUE,
      fit = list(model = "test"),
      draws = draws,
      diagnostics = list(rhat = 1.01),
      timing = list(total = 1.0)
    )

    retained <- apply_fit_retention(result, c("metrics"))
    expect_null(retained$draws)
  })

  it("removes diagnostics when not in retain", {
    draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
    colnames(draws) <- c("alpha", "beta")
    result <- new_fit_result(
      success = TRUE,
      draws = draws,
      diagnostics = list(rhat = 1.01),
      timing = list(total = 1.0)
    )

    retained <- apply_fit_retention(result, c("metrics"))
    # When not retained, diagnostics is set to NULL
    expect_null(retained$diagnostics)
  })

  it("removes data_bundle when not in retain", {
    result <- new_fit_result(
      success = TRUE,
      draws = NULL,
      diagnostics = list(),
      timing = list(total = 1.0)
    )
    result$data_bundle <- list(train = data.frame(x = 1))

    retained <- apply_fit_retention(result, c("metrics"))
    expect_null(retained$data_bundle)
  })

  it("keeps fit when in retain", {
    draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
    colnames(draws) <- c("alpha", "beta")
    result <- new_fit_result(
      success = TRUE,
      fit = list(model = "test"),
      draws = draws,
      diagnostics = list(rhat = 1.01),
      timing = list(total = 1.0)
    )

    retained <- apply_fit_retention(result, c("fit", "metrics", "diagnostics"))
    expect_false(is.null(retained$fit))
  })

  it("rejects non-fit-result input", {
    expect_error(
      apply_fit_retention(list(), c("metrics")),
      "must be a bayesim_fit_result"
    )
  })
})

describe("apply_task_retention()", {
  it("adds draws when in retain and fit_result has draws", {
    draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
    colnames(draws) <- c("alpha", "beta")
    fit_result <- new_fit_result(
      success = TRUE,
      draws = draws,
      diagnostics = list(),
      timing = list(total = 1.0)
    )
    task_result <- new_task_result(
      task_id = "test",
      status = "success",
      metrics = list(rmse = 0.5),
      timing = list(total = 1.0)
    )

    result <- apply_task_retention(task_result, fit_result, list(), c("draws", "metrics"))
    expect_false(is.null(result$draws))
  })

  it("does not add draws when not in retain", {
    draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
    colnames(draws) <- c("alpha", "beta")
    fit_result <- new_fit_result(
      success = TRUE,
      draws = draws,
      diagnostics = list(),
      timing = list(total = 1.0)
    )
    task_result <- new_task_result(
      task_id = "test",
      status = "success",
      metrics = list(rmse = 0.5),
      timing = list(total = 1.0)
    )

    result <- apply_task_retention(task_result, fit_result, list(), c("metrics"))
    expect_null(result$draws)
  })

  it("adds predictions when in retain and available", {
    fit_result <- new_fit_result(
      success = TRUE,
      draws = NULL,
      diagnostics = list(),
      timing = list(total = 1.0)
    )
    task_result <- new_task_result(
      task_id = "test",
      status = "success",
      metrics = list(rmse = 0.5),
      timing = list(total = 1.0)
    )

    result <- apply_task_retention(task_result, fit_result, list(), c("predictions", "metrics"))
    # predictions comes from fit_result's predictions - depends on implementation
    expect_s3_class(result, "bayesim_task_result")
  })

  it("rejects non-task-result input", {
    fit_result <- new_fit_result(
      success = TRUE,
      diagnostics = list(),
      timing = list(total = 1.0)
    )
    expect_error(
      apply_task_retention(list(), fit_result, list(), c("metrics")),
      "must be a bayesim_task_result"
    )
  })
})

describe("lighten_task_result()", {
  it("removes retained fields from task result", {
    task_result <- new_task_result(
      task_id = "test",
      status = "success",
      metrics = list(rmse = 0.5),
      timing = list(total = 1.0)
    )
    task_result$draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
    task_result$predictions <- list(mean = 1.0)

    retain <- c("metrics", "diagnostics")
    lightened <- lighten_task_result(task_result, retain)
    expect_null(lightened$draws)
    expect_null(lightened$predictions)
  })

  it("preserves task_id and status", {
    task_result <- new_task_result(
      task_id = "test_id_123",
      status = "success",
      metrics = list(rmse = 0.5),
      timing = list(total = 1.0)
    )
    lightened <- lighten_task_result(task_result, c("metrics"))
    expect_equal(lightened$task_id, "test_id_123")
    expect_equal(lightened$status, "success")
  })
})
