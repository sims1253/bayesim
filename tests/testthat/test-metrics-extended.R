# Tests for the extended metric library (M5). Uses MockFitter-style fixtures
# rather than real brms compiles so the tests are fast and CRAN-safe.

# Build a synthetic fit_result + data_bundle + context for metric.compute().
make_fixture <- function(draws = NULL,
                         true_params = NULL,
                         vars_of_interest = NULL,
                         diagnostics = NULL,
                         predictions = NULL,
                         log_lik = NULL,
                         loo = NULL,
                         test = NULL,
                         response = "y") {
  if (is.null(draws)) {
    draws <- matrix(c(rnorm(200, mean = 1), rnorm(200, mean = 2)), ncol = 2,
      dimnames = list(NULL, c("b_x", "b_Intercept")))
  }
  if (is.null(true_params)) true_params <- c(b_x = 1, b_Intercept = 2)
  if (is.null(vars_of_interest)) vars_of_interest <- c("b_x", "b_Intercept")
  list(
    fit_result = list(draws = draws, diagnostics = diagnostics),
    data_bundle = list(
      train = data.frame(y = rnorm(10), x = rnorm(10)),
      test = test,
      response = response,
      true_params = true_params,
      vars_of_interest = vars_of_interest
    ),
    context = list(predictions = predictions, log_lik = log_lik, loo = loo),
    task_ctx = list(task_id = "t1", rep_idx = 1L)
  )
}

describe("M5 extended metrics", {
  it("pred_mae_metric computes mean absolute error", {
    fx <- make_fixture(predictions = list(predicted_mean = 1:10))
    fx$data_bundle$response <- "y"
    fx$data_bundle$train$y <- 2:11
    out <- compute_metric(pred_mae_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_equal(out$value, 1)
    expect_equal(out$n_obs, 10)
  })

  it("pred_mae_metric returns NA when predictions absent", {
    fx <- make_fixture(predictions = NULL)
    out <- compute_metric(pred_mae_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_true(is.na(out$value))
  })

  it("pred_mse_metric computes mean squared error", {
    fx <- make_fixture(predictions = list(predicted_mean = rep(0, 5)))
    fx$data_bundle$train$y <- c(1, 2, 3, 4, 5)
    out <- compute_metric(pred_mse_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_equal(out$value, mean((1:5)^2))
  })

  it("pos_prob_metric computes fraction of positive draws", {
    draws <- matrix(c(rep(1, 80), rep(-1, 120)), ncol = 1, dimnames = list(NULL, "b_x"))
    fx <- make_fixture(draws = draws, true_params = c(b_x = 0), vars_of_interest = "b_x")
    out <- compute_metric(pos_prob_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_equal(out$by_param[["b_x"]], 0.4)
    expect_equal(out$mean, 0.4)
  })

  it("posterior_summary_metric returns mean/median/sd/quantiles", {
    fx <- make_fixture()
    out <- compute_metric(posterior_summary_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_true(all(c("b_x", "b_Intercept") %in% names(out$mean)))
    expect_equal(out$mean[["b_x"]], mean(fx$fit_result$draws[, "b_x"]))
    expect_equal(out$median[["b_Intercept"]], median(fx$fit_result$draws[, "b_Intercept"]))
    expect_true(out$q_lower[["b_x"]] < out$q_upper[["b_x"]])
  })

  it("convergence_metric surfaces diagnostics", {
    fx <- make_fixture(diagnostics = list(
      rhat_max = 1.01, ess_bulk_min = 100, ess_tail_min = 90, divergent = 0L
    ))
    out <- compute_metric(convergence_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_equal(out$rhat_max, 1.01)
    expect_equal(out$ess_bulk_min, 100)
    expect_equal(out$divergent, 0L)
  })

  it("convergence_metric returns NAs when diagnostics absent", {
    fx <- make_fixture(diagnostics = NULL)
    out <- compute_metric(convergence_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_true(is.na(out$rhat_max))
  })

  it("sampler_diagnostics_metric surfaces divergences/treedepth", {
    fx <- make_fixture(diagnostics = list(divergent = 3L, max_treedepth = 5L))
    out <- compute_metric(sampler_diagnostics_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_equal(out$divergent, 3L)
    expect_equal(out$max_treedepth, 5L)
  })

  it("rank_metric computes SBC ranks (draws < true value)", {
    # 200 draws of b_x centered at 1; true b_x = 0 -> all draws above -> rank 0
    draws <- matrix(rnorm(200, mean = 5, sd = 0.1), ncol = 1, dimnames = list(NULL, "b_x"))
    fx <- make_fixture(draws = draws, true_params = c(b_x = 0), vars_of_interest = "b_x")
    out <- compute_metric(rank_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_equal(out$by_param[["b_x"]], 0L)
    expect_equal(out$n_draws, 200)
    # F4: n_ranks per variable present (post-thinning sample size + 1).
    expect_true("b_x" %in% names(out$n_ranks))
    expect_true(out$n_ranks[["b_x"]] >= 2L)
  })

  it("rank_metric thin=FALSE disables thinning (stride 1)", {
    draws <- matrix(rnorm(200, mean = 5, sd = 0.1), ncol = 1, dimnames = list(NULL, "b_x"))
    fx <- make_fixture(draws = draws, true_params = c(b_x = 0), vars_of_interest = "b_x")
    out <- compute_metric(rank_metric(thin = FALSE), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_equal(out$stride, 1L)
    expect_equal(out$n_ranks[["b_x"]], 201L)  # 200 draws + 1 possible ranks
  })

  it("rank_metric integer thin uses the stride directly", {
    draws <- matrix(rnorm(300, mean = 5, sd = 0.1), ncol = 1, dimnames = list(NULL, "b_x"))
    fx <- make_fixture(draws = draws, true_params = c(b_x = 0), vars_of_interest = "b_x")
    out <- compute_metric(rank_metric(thin = 10L), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_equal(out$stride, 10L)
    expect_equal(out$n_ranks[["b_x"]], 31L)  # floor(300/10)=30 kept + 1
  })

  it("rank_metric auto-thinning reduces n_ranks on autocorrelated draws", {
    # AR(1) draws: high autocorrelation -> low ESS -> auto-thinning keeps fewer.
    set.seed(42)
    n <- 1000
    ar <- 0.95
    x <- numeric(n); x[1] <- rnorm(1)
    for (i in seq_len(n - 1L)) x[i + 1L] <- ar * x[i] + rnorm(1)
    draws <- matrix(x, ncol = 1, dimnames = list(NULL, "b_x"))
    fx <- make_fixture(draws = draws, true_params = c(b_x = 0), vars_of_interest = "b_x")
    out_auto <- compute_metric(rank_metric(thin = "auto"), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    out_none <- compute_metric(rank_metric(thin = FALSE), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_true(out_auto$stride > 1L)
    expect_true(out_auto$n_ranks[["b_x"]] < out_none$n_ranks[["b_x"]])
  })

  it("rank_metric output passes flatten with 2 params (named vectors, F4)", {
    # The flatten/validate layer requires multi-element named vectors to be
    # double AND named; verify a 2-param rank result flattens to the expected
    # column names without error.
    draws <- matrix(rnorm(400), ncol = 2, dimnames = list(NULL, c("b_x", "b_Intercept")))
    fx <- make_fixture(draws = draws,
                       true_params = c(x = 0, Intercept = 0),
                       vars_of_interest = c("x", "Intercept"))
    out <- compute_metric(rank_metric(thin = FALSE), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    flat <- flatten_metric_output(out, "rank")
    expect_true("rank__by_param__x" %in% names(flat))
    expect_true("rank__by_param__Intercept" %in% names(flat))
    expect_true("rank__n_ranks__x" %in% names(flat))
    expect_true("rank__n_ranks__Intercept" %in% names(flat))
  })

  it("rank_metric ranks are ~uniform on 0..S for independent draws (chi-square)", {
    # Under correct calibration with independent draws, ranks of a truth drawn
    # from the same distribution are uniform on 0..S. Use a generous alpha.
    skip_if_not(requireNamespace("posterior", quietly = TRUE))
    set.seed(123)
    S <- 200L
    n_reps <- 200L
    truth <- rnorm(n_reps)
    ranks <- vapply(truth, function(tv) {
      draws <- matrix(rnorm(S), ncol = 1, dimnames = list(NULL, "b_x"))
      fx <- make_fixture(draws = draws, true_params = c(b_x = tv), vars_of_interest = "b_x")
      out <- compute_metric(rank_metric(thin = FALSE), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
      as.integer(out$by_param[["b_x"]])
    }, integer(1))
    # Chi-square uniformity: bin into 10 bins over 0..S.
    bins <- cut(ranks, breaks = seq(0, S + 1L, length.out = 11L), include.lowest = TRUE)
    obs <- as.numeric(table(bins))
    exp <- rep(n_reps / 10, 10)
    chisq <- sum((obs - exp)^2 / exp)
    # df=9; generous alpha 0.001 critical value ~ 27.88.
    expect_true(chisq < 30,
      info = paste("chi-square", round(chisq, 2), "exceeds uniformity threshold"))
  })

  it("rank_metric returns no ranks when true_params absent (E5: no mean field)", {
    # Build a fixture with explicitly NULL true_params (bypass make_fixture defaults).
    draws <- matrix(rnorm(200), ncol = 1, dimnames = list(NULL, "b_x"))
    fit_result <- list(draws = draws, diagnostics = NULL)
    data_bundle <- list(
      train = data.frame(y = rnorm(10), x = rnorm(10)), test = NULL,
      response = "y", true_params = NULL, vars_of_interest = "b_x"
    )
    out <- compute_metric(rank_metric(), fit_result, data_bundle, list(), list(task_id = "t1"))
    expect_null(out$mean)
    expect_true(is.na(out$n_draws))
  })

  it("elpd_loo_metric surfaces loo context fields", {
    fx <- make_fixture(loo = list(elpd = -10.5, p_loo = 2.1, elpd_se = 1.3, pareto_k = c(0.1, 0.7)))
    out <- compute_metric(elpd_loo_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_equal(out$elpd, -10.5)
    expect_equal(out$p_loo, 2.1)
    expect_equal(out$pareto_k_max, 0.7)
  })

  it("elpd_loo_metric returns NAs when loo absent", {
    fx <- make_fixture(loo = NULL)
    out <- compute_metric(elpd_loo_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_true(is.na(out$elpd))
  })

  it("rmse_loo_metric degrades to NA without PSIS/epred context (F3)", {
    # The legacy elpd-only proxy formula was removed in F3; rmse_loo now needs
    # the PSIS object + pointwise log_lik + a prediction matrix from the
    # metric context. With only the loo summary present, it must NA gracefully.
    fx <- make_fixture(loo = list(elpd = -6, pointwise = cbind(elpd_loo = c(-1, -2, -3))))
    out <- compute_metric(rmse_loo_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_true(is.na(out$value))
    expect_true(is.na(out$pareto_k_max))
  })

  it("r2_loo_metric degrades to NA without PSIS/epred context (F3)", {
    fx <- make_fixture(loo = list(elpd = -10, pointwise = cbind(elpd_loo = rep(-1, 10))))
    out <- compute_metric(r2_loo_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_true(is.na(out$value))
  })

  it("elpd_test_metric computes log-sum-exp elpd on log_lik", {
    ll <- matrix(c(-1, -2, -1.5, -2.5), nrow = 2) # 2 draws x 2 obs
    fx <- make_fixture(log_lik = ll, test = data.frame(y = c(1, 2)))
    out <- compute_metric(elpd_test_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    # obs1: max(-1,-1.5)=-1, log(mean(exp(0,-.5)))=log(mean(1,0.607))=log(0.803)=-0.219
    expect_equal(out$n_obs, 2)
    expect_true(is.numeric(out$value))
  })

  it("elpd_test_metric returns NA when test absent", {
    fx <- make_fixture(test = NULL)
    out <- compute_metric(elpd_test_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_true(is.na(out$value))
  })

  it("rmse_test_metric computes RMSE on test set", {
    fx <- make_fixture(
      predictions = list(predicted_mean = c(1, 2, 3)),
      test = data.frame(y = c(2, 2, 2))
    )
    out <- compute_metric(rmse_test_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_equal(out$value, sqrt(mean((c(1, 2, 3) - c(2, 2, 2))^2)))
    expect_equal(out$n_obs, 3)
  })

  it("r2_test_metric computes R-squared on test set", {
    fx <- make_fixture(
      predictions = list(predicted_mean = c(1, 2, 3)),
      test = data.frame(y = c(1, 2, 3))
    )
    out <- compute_metric(r2_test_metric(), fx$fit_result, fx$data_bundle, fx$context, fx$task_ctx)
    expect_equal(out$value, 1) # perfect prediction
  })
})
