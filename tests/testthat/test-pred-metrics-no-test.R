# E2: prediction metrics refuse to silently fall back to the training set.
skip_on_cran()

describe("E2: prediction metrics require a test set", {
  it("pred_rmse_metric returns NA when no test set", {
    draws <- matrix(
      rnorm(200),
      ncol = 2,
      dimnames = list(NULL, c("b_x", "sigma"))
    )
    fit_result <- list(draws = draws)
    data_bundle <- list(
      train = data.frame(y = rnorm(10), x = rnorm(10)),
      test = NULL,
      response = "y"
    )
    context <- list(predictions = list(predicted_mean = rnorm(10)))
    out <- compute_metric(
      pred_rmse_metric(),
      fit_result,
      data_bundle,
      context,
      list(task_id = "t")
    )
    expect_true(is.na(out$value))
    expect_true(is.na(out$n_obs))
  })

  it("pred_mae_metric computes on the test set when present", {
    actual <- c(1, 2, 3, 4)
    data_bundle <- list(
      train = data.frame(y = 0, x = 0),
      test = data.frame(y = actual, x = c(0, 0, 0, 0)),
      response = "y"
    )
    fit_result <- list(draws = matrix(rnorm(20), ncol = 2))
    context <- list(predictions = list(predicted_mean = c(1.1, 1.9, 3.2, 3.8)))
    out <- compute_metric(
      pred_mae_metric(),
      fit_result,
      data_bundle,
      context,
      list(task_id = "t")
    )
    expect_false(is.na(out$value))
    expect_equal(out$n_obs, 4L)
  })
})
