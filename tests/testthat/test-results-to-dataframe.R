# test-results-to-dataframe.R
# Tests for results_to_dataframe conversion logic

library(bayesim)

describe("results_to_dataframe()", {
  it("returns empty data frame for NULL input", {
    result <- results_to_dataframe(NULL)
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 0)
  })

  it("returns empty data frame for empty list", {
    result <- results_to_dataframe(list())
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 0)
  })

  it("rejects non-list input", {
    expect_error(
      results_to_dataframe("not a list"),
      "must be a list"
    )
  })

  it("converts single successful task result", {
    task_result <- new_task_result(
      task_id = "d001_f001_r00001",
      status = "success",
      metrics = list(rmse = 0.5, bias = 0.01),
      diagnostics = list(n_eff = 500, rhat = 1.01),
      timing = list(total = 5.2)
    )

    result <- results_to_dataframe(list(task_result))
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 1)
    expect_equal(result$task_id, "d001_f001_r00001")
    expect_equal(result$status, "success")
    expect_equal(result$rmse, 0.5)
    expect_equal(result$bias, 0.01)
  })

  it("converts single failed task result", {
    task_result <- new_task_result(
      task_id = "d001_f001_r00002",
      status = "failed",
      error = list(
        error_class = "convergence_error",
        error_message = "R-hat > 1.1"
      ),
      timing = list(total = 2.0)
    )

    result <- results_to_dataframe(list(task_result))
    expect_equal(result$error_class, "convergence_error")
    expect_equal(result$error_message, "R-hat > 1.1")
  })

  it("skips NULL entries", {
    task_result <- new_task_result(
      task_id = "d001_f001_r00001",
      status = "success",
      metrics = list(rmse = 0.5),
      timing = list(total = 1.0)
    )

    result <- results_to_dataframe(list(NULL, task_result, NULL))
    expect_equal(nrow(result), 1)
  })

  it("handles multiple task results with different metric columns", {
    t1 <- new_task_result(
      task_id = "d001_f001_r00001",
      status = "success",
      metrics = list(rmse = 0.5),
      timing = list(total = 1.0)
    )
    t2 <- new_task_result(
      task_id = "d001_f001_r00002",
      status = "success",
      metrics = list(rmse = 0.3, bias = 0.01),
      timing = list(total = 2.0)
    )

    result <- results_to_dataframe(list(t1, t2))
    expect_equal(nrow(result), 2)
    # t2 has bias but t1 doesn't - should be NA for t1
    expect_true(is.na(result$bias[1]))
    expect_equal(result$bias[2], 0.01)
  })

  it("flattens named diagnostic vectors", {
    task_result <- new_task_result(
      task_id = "test",
      status = "success",
      metrics = list(rmse = 0.5),
      diagnostics = list(rhat = c(alpha = 1.01, beta = 1.00)),
      timing = list(total = 1.0)
    )

    result <- results_to_dataframe(list(task_result))
    expect_true("rhat__alpha" %in% names(result))
    expect_true("rhat__beta" %in% names(result))
  })

  it("includes timing_total when available", {
    task_result <- new_task_result(
      task_id = "test",
      status = "success",
      metrics = list(rmse = 0.5),
      timing = list(total = 3.14)
    )

    result <- results_to_dataframe(list(task_result))
    expect_equal(result$timing_total, 3.14)
  })
})
