# D1: LinearRegressionFitter — conjugate NIG Bayesian linear regression.
# Tests: analytic posterior mean/cov vs closed form; coverage ~0.95 over
# replicates; validate_fitter smoke test; full contract orientation (S x N).

skip_on_cran()

.gen <- function(data_spec, task_ctx) {
  n <- data_spec$n %||% 200L
  beta <- data_spec$beta %||% 0.5
  x <- stats::rnorm(n)
  y <- 2 + beta * x + stats::rnorm(n)
  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(Intercept = 2, x = beta, sigma = 1),
    vars_of_interest = c("Intercept", "x", "sigma"),
    references = c(Intercept = 0, x = 0, sigma = 1),
    meta = list()
  )
}

describe("LinearRegressionFitter analytic posterior", {
  it("posterior means match OLS (weak prior)", {
    withr::local_seed(1)
    n <- 400L
    x <- stats::rnorm(n)
    y <- 1 + 0.7 * x + stats::rnorm(n, sd = 0.5)
    data_bundle <- list(
      train = data.frame(y = y, x = x),
      test = NULL,
      response = "y",
      true_params = c(Intercept = 1, x = 0.7, sigma = 0.5),
      vars_of_interest = c("Intercept", "x", "sigma")
    )
    fitter <- LinearRegressionFitter(n_draws = 8000L)
    res <- fit_model(fitter, data_bundle, list(formula = y ~ x), seed = 42L,
                     task_ctx = list(task_id = "t"))
    draws <- res$draws
    ols <- lm(y ~ x)
    expect_equal(colnames(draws), c("Intercept", "x", "sigma"))
    expect_equal(mean(draws[, "Intercept"]), unname(coef(ols)[1]), tolerance = 0.05)
    expect_equal(mean(draws[, "x"]), unname(coef(ols)[2]), tolerance = 0.05)
    expect_equal(mean(draws[, "sigma"]), summary(ols)$sigma, tolerance = 0.05)
  })

  it("posterior covariance matches closed-form NIG", {
    # With a known weak prior, the posterior covariance of beta is
    # sigma^2 * (X'X)^{-1}; check the (x,x) entry against the draws.
    withr::local_seed(7)
    n <- 300L
    x <- stats::rnorm(n)
    X <- cbind(1, x)
    y <- as.vector(X %*% c(0, 1) + stats::rnorm(n))
    data_bundle <- list(
      train = data.frame(y = y, x = x), test = NULL, response = "y",
      true_params = c(Intercept = 0, x = 1, sigma = 1),
      vars_of_interest = c("Intercept", "x", "sigma")
    )
    fitter <- LinearRegressionFitter(n_draws = 20000L)
    res <- fit_model(fitter, data_bundle, list(formula = y ~ x), seed = 9L,
                     task_ctx = list(task_id = "t"))
    sigma2_hat <- mean(res$draws[, "sigma"]^2)
    closed_var_x <- sigma2_hat * solve(crossprod(X))[2, 2]
    expect_equal(stats::var(res$draws[, "x"]), closed_var_x, tolerance = 0.1)
  })
})

describe("LinearRegressionFitter contract", {
  it("passes validate_fitter smoke test", {
    expect_silent(
      validate_fitter(LinearRegressionFitter(n_draws = 100L),
                      smoke_test = TRUE, verbose = FALSE)
    )
  })

  it("log_lik / epred / predicted_samples are S x N", {
    withr::local_seed(3)
    n <- 25L
    data_bundle <- list(
      train = data.frame(y = stats::rnorm(n), x = stats::rnorm(n)),
      test = NULL, response = "y"
    )
    fitter <- LinearRegressionFitter(n_draws = 100L)
    res <- fit_model(fitter, data_bundle, list(formula = y ~ x), seed = 1L,
                     task_ctx = list(task_id = "t"))

    ll <- log_lik_matrix(fitter, res)
    expect_equal(dim(ll), c(100L, n))

    epred <- predict_epred(fitter, res)
    expect_equal(dim(epred), c(100L, n))

    preds <- predict_fit(fitter, res)
    expect_equal(dim(preds$predicted_samples), c(100L, n))
    expect_equal(length(preds$predicted_mean), n)
    expect_equal(length(preds$predicted_sd), n)
  })

  it("fit_diagnostics reports rhat 1, ESS = n_draws, 0 divergences", {
    data_bundle <- list(
      train = data.frame(y = 1:10, x = 1:10), test = NULL, response = "y"
    )
    fitter <- LinearRegressionFitter(n_draws = 250L)
    res <- fit_model(fitter, data_bundle, list(formula = y ~ x), seed = 1L,
                     task_ctx = list(task_id = "t"))
    diag <- fit_diagnostics(fitter, res)
    expect_equal(diag$rhat_max, 1)
    expect_equal(diag$ess_bulk, 250L)
    expect_equal(diag$divergent, 0L)
  })
})

describe("LinearRegressionFitter coverage", {
  it("coverage ~ 0.95 over replicates", {
    skip_on_cran()
    # Per-task 90% intervals of 'x'; coverage should be near 0.90. We use a
    # generous band given MC noise on a small number of replicates.
    fitter <- LinearRegressionFitter(n_draws = 500L)
    n_rep <- 200L
    config <- simulation_config(
      data_grid = data.frame(n = 100L, beta = 0.5),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = fitter,
      metrics = list(coverage_metric()),
      n_replicates = n_rep,
      seed = 2024L
    )
    # coverage_metric needs posterior_summary-style intervals; supply both.
    config@metrics <- list(coverage_metric(), posterior_summary_metric())

    result <- run_simulation(config, resume = "never", progress = FALSE)
    # All tasks succeeded.
    expect_true(all(result$summary$status == "success"))

    # Compute coverage of the 'x' coefficient from posterior summaries:
    # count replicates whose truth (beta=0.5) lies in the 90% posterior interval.
    lo <- result$summary[["posterior_summary__q_lower__x"]]
    hi <- result$summary[["posterior_summary__q_upper__x"]]
    covered <- !is.na(lo) & !is.na(hi) & (0.5 >= lo) & (0.5 <= hi)
    cov <- mean(covered)
    # Nominal 0.90; allow a generous Monte-Carlo band.
    expect_true(cov > 0.80 && cov < 0.97)
  })
})
