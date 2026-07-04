# test-fitter-orientation.R
# Tests for the S x N (draws x observations) log_lik matrix convention (step A1).
#
# The brms/loo convention is that log_lik() returns an S x N matrix: S rows
# (one per posterior draw) and N columns (one per observation). These tests
# verify (a) MockFitter honors the convention, (b) validate_fitter() rejects a
# fitter that returns the transposed (N x S) orientation, and (c) downstream
# consumers (elpd_test_metric) read N from the columns correctly.

library(bayesim)

describe("log_lik S x N orientation convention", {
  it("MockFitter log_lik returns S x N (draws x observations)", {
    fitter <- MockFitter()
    n_obs <- 15L
    data_bundle <- list(
      train = data.frame(x = seq_len(n_obs), y = rnorm(n_obs)),
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
    # N observations => N columns.
    expect_equal(ncol(ll), n_obs)
    # S draws => S rows. Default MockFitter: n_draws=100, n_chains=4 => 400.
    s_expected <- as.integer(fitter@n_draws) * as.integer(fitter@n_chains)
    expect_equal(nrow(ll), s_expected)
    # Sanity: n_obs and n_draws differ so orientation is unambiguous.
    expect_true(s_expected != n_obs)
  })

  it("validate_fitter passes for the correctly-oriented MockFitter", {
    # The default MockFitter honors S x N, so the smoke test should pass.
    expect_silent(
      validate_fitter(MockFitter(), smoke_test = TRUE, verbose = FALSE)
    )
  })

  it("validate_fitter rejects a fitter whose log_lik is transposed (N x S)", {
    # A custom Fitter subclass that implements every method the smoke test
    # exercises, but deliberately returns log_lik as N x S (observations in
    # rows, draws in columns) -- the WRONG orientation.
    TransposedFitter <- S7::new_class(
      "TransposedFitter",
      parent = Fitter,
      properties = list(
        name = S7::new_property(S7::class_character, default = "transposed"),
        supports_predictions = S7::new_property(S7::class_logical, default = TRUE),
        supports_log_lik = S7::new_property(S7::class_logical, default = TRUE),
        supports_loo = S7::new_property(S7::class_logical, default = TRUE),
        n_draws = S7::new_property(S7::class_integer, default = 50L)
      )
    )

    S7::method(fit, TransposedFitter) <- function(
      fitter, data_bundle, fit_spec, seed, task_ctx
    ) {
      n_obs <- nrow(data_bundle$train)
      n_draws <- as.integer(fitter@n_draws)
      draws <- matrix(rnorm(n_draws * 3), nrow = n_draws, ncol = 3)
      colnames(draws) <- c("intercept", "beta", "sigma")
      mock_fit <- list(
        fitter = "transposed",
        data_bundle = data_bundle,
        seed = seed %||% 12345L,
        n_obs = n_obs
      )
      new_fit_result(
        success = TRUE,
        fit = mock_fit,
        draws = draws,
        diagnostics = list(rhat_max = 1.0, ess_bulk = 100, ess_tail = 100, divergent = 0L),
        timing = list(total = 1.0)
      )
    }

    S7::method(extract_draws, TransposedFitter) <- function(fitter, fit_result, variables = NULL) {
      draws <- fit_result$draws
      if (!is.null(variables)) draws <- draws[, variables, drop = FALSE]
      draws
    }

    S7::method(predict_fit, TransposedFitter) <- function(
      fitter, fit_result, newdata = NULL, seed = NULL
    ) {
      data_bundle <- fit_result$fit$data_bundle
      data <- newdata %||% data_bundle$train
      n_obs <- nrow(data)
      n_draws <- as.integer(fitter@n_draws)
      predicted_samples <- matrix(rnorm(n_obs * n_draws), nrow = n_obs, ncol = n_draws)
      list(
        predicted_mean = rowMeans(predicted_samples),
        predicted_samples = predicted_samples,
        predicted_sd = apply(predicted_samples, 1, sd)
      )
    }

    S7::method(log_lik, TransposedFitter) <- function(fitter, fit_result, newdata = NULL) {
      # DELIBERATELY WRONG: N x S (observations in rows, draws in columns).
      data_bundle <- fit_result$fit$data_bundle
      data <- newdata %||% data_bundle$train
      n_obs <- nrow(data)
      n_draws <- as.integer(fitter@n_draws)
      matrix(rnorm(n_obs * n_draws, mean = -2, sd = 0.5), nrow = n_obs, ncol = n_draws)
    }

    S7::method(loo_fit, TransposedFitter) <- function(fitter, fit_result) {
      list(elpd = -10, p_loo = 1, elpd_se = 1, pareto_k = numeric())
    }

    S7::method(diagnostics, TransposedFitter) <- function(fitter, fit_result) {
      list(rhat_max = 1.0, ess_bulk = 100, ess_tail = 100, divergent = 0L)
    }

    # n_obs in the smoke test is 20; n_draws here is 50, so ncol(ll) = 50 != 20
    # must trigger the orientation check.
    expect_error(
      validate_fitter(TransposedFitter(), smoke_test = TRUE, verbose = FALSE),
      class = "bayesim_validation_error"
    )
    expect_error(
      validate_fitter(TransposedFitter(), smoke_test = TRUE, verbose = FALSE),
      "wrong number of columns"
    )
  })

  it("elpd_test_metric reports n_obs from log_lik columns (S x N)", {
    # Build a real fit_result with a non-square (n_obs != n_draws) log_lik via
    # MockFitter, then feed it to elpd_test_metric through compute().
    fitter <- MockFitter()
    n_obs <- 18L
    data_bundle <- list(
      train = data.frame(x = seq_len(n_obs), y = rnorm(n_obs)),
      test = data.frame(x = seq_len(n_obs), y = rnorm(n_obs)),
      response = "y"
    )
    fit_result <- fit(
      fitter,
      data_bundle,
      fit_spec = data.frame(model = "test"),
      seed = 7L,
      task_ctx = list(task_id = "test")
    )

    ll <- log_lik(fitter, fit_result)
    # Orientation sanity before computing the metric.
    expect_equal(ncol(ll), n_obs)
    s_expected <- as.integer(fitter@n_draws) * as.integer(fitter@n_chains)
    expect_true(nrow(ll) != n_obs)

    context <- list(log_lik = ll)
    out <- compute(
      elpd_test_metric(),
      fit_result,
      data_bundle,
      context,
      list(task_id = "t1")
    )
    # n_obs is read from ncol(ll) under the S x N convention.
    expect_equal(out$n_obs, n_obs)
    expect_true(is.numeric(out$value))
    expect_true(is.finite(out$value))
  })
})
