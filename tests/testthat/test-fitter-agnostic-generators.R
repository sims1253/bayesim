# Workstream I1: fitter-agnostic generators (prior_draws_generator,
# forward_sim_generator). Tested with LinearRegressionFitter (no Stan).

# A predictor generator that consumes the ambient RNG state.
gaussian_predictors <- function(data_spec, task_ctx) {
  n <- as.integer(data_spec$n %||% 20L)
  data.frame(x = stats::rnorm(n))
}

# A small pilot data_bundle used for the one-time preconditioning fit.
pilot_bundle <- list(
  train = data.frame(y = stats::rnorm(30), x = stats::rnorm(30)),
  test = NULL,
  response = "y"
)

describe("prior_draws_generator", {
  it("returns a valid data_bundle for LinearRegressionFitter", {
    withr::local_seed(42)
    fitter <- LinearRegressionFitter(n_draws = 50L)
    gen <- prior_draws_generator(
      fitter = fitter,
      fit_spec = list(formula = y ~ x),
      pilot_bundle = pilot_bundle,
      predictor_generator = gaussian_predictors
    )

    bundle <- gen(list(n = 20L), list(task_id = "t1", rep_idx = 1L))

    # Must pass validate_data_bundle.
    expect_silent(bayesim:::validate_data_bundle(bundle))

    # Structure.
    expect_s3_class(bundle$train, "data.frame")
    expect_equal(bundle$response, "y")
    expect_null(bundle$test)
    expect_equal(bundle$meta$generator, "prior_draws")
    expect_equal(bundle$meta$truth_draw_id, 1L)
    # Non-brms path degrades to posterior draws (approximate prior).
    expect_equal(bundle$meta$prior_source, "posterior_degraded")
  })

  it("true_params names match the fitter's draw column names", {
    withr::local_seed(1)
    fitter <- LinearRegressionFitter(n_draws = 30L)
    gen <- prior_draws_generator(
      fitter = fitter,
      fit_spec = list(formula = y ~ x),
      pilot_bundle = pilot_bundle,
      predictor_generator = gaussian_predictors
    )
    bundle <- gen(list(n = 20L), list(task_id = "t1", rep_idx = 1L))

    # LinearRegressionFitter draw columns: Intercept, x, sigma.
    expect_setequal(names(bundle$true_params), c("Intercept", "x", "sigma"))
    expect_setequal(bundle$vars_of_interest, c("Intercept", "x", "sigma"))
    # names(true_params) must match vars_of_interest exactly.
    expect_equal(names(bundle$true_params), bundle$vars_of_interest)
  })

  it("pins theta deterministically by rep_idx and simulates a fresh y", {
    withr::local_seed(2)
    fitter <- LinearRegressionFitter(n_draws = 40L)
    gen <- prior_draws_generator(
      fitter = fitter,
      fit_spec = list(formula = y ~ x),
      pilot_bundle = pilot_bundle,
      predictor_generator = gaussian_predictors
    )
    b1 <- gen(list(n = 20L), list(task_id = "a", rep_idx = 1L))
    b2 <- gen(list(n = 20L), list(task_id = "b", rep_idx = 1L))
    expect_equal(b1$true_params, b2$true_params)

    # Different rep_idx -> different theta.
    b3 <- gen(list(n = 20L), list(task_id = "c", rep_idx = 2L))
    expect_false(identical(b1$true_params, b3$true_params))
    expect_equal(b3$meta$truth_draw_id, 2L)
  })

  it("wraps rep_idx modulo n_draws without indexing errors", {
    withr::local_seed(3)
    n_draws <- 25L
    fitter <- LinearRegressionFitter(n_draws = n_draws)
    gen <- prior_draws_generator(
      fitter = fitter,
      fit_spec = list(formula = y ~ x),
      pilot_bundle = pilot_bundle,
      predictor_generator = gaussian_predictors
    )
    rep_idx_wrap <- n_draws + 5L
    expected_id <- ((rep_idx_wrap - 1L) %% n_draws) + 1L
    bw <- gen(list(n = 20L), list(task_id = "tw", rep_idx = rep_idx_wrap))
    expect_equal(bw$meta$truth_draw_id, expected_id)
  })

  it("rejects a non-Fitter object", {
    expect_error(
      prior_draws_generator(
        fitter = "not_a_fitter",
        fit_spec = list(formula = y ~ x),
        pilot_bundle = pilot_bundle,
        predictor_generator = gaussian_predictors
      ),
      class = "bayesim_config_error"
    )
  })

  it("rejects a missing pilot_bundle$train", {
    expect_error(
      prior_draws_generator(
        fitter = LinearRegressionFitter(n_draws = 10L),
        fit_spec = list(formula = y ~ x),
        pilot_bundle = list(),
        predictor_generator = gaussian_predictors
      ),
      class = "bayesim_config_error"
    )
  })
})

describe("forward_sim_generator", {
  it("returns a valid data_bundle for LinearRegressionFitter", {
    withr::local_seed(11)
    fitter <- LinearRegressionFitter(n_draws = 50L)
    gen <- forward_sim_generator(
      fitter = fitter,
      fit_spec = list(formula = y ~ x),
      pilot_bundle = pilot_bundle,
      predictor_generator = gaussian_predictors
    )

    bundle <- gen(list(n = 20L), list(task_id = "t1", rep_idx = 1L))
    expect_silent(bayesim:::validate_data_bundle(bundle))

    expect_s3_class(bundle$train, "data.frame")
    expect_true("y" %in% names(bundle$train))
    expect_false(all(is.na(bundle$train$y)))
    expect_equal(bundle$response, "y")
    expect_null(bundle$test)
    expect_equal(bundle$meta$generator, "forward_sim")
    expect_equal(bundle$meta$truth_draw_id, 1L)
  })

  it("true_params names match the fitter's draw column names", {
    withr::local_seed(12)
    fitter <- LinearRegressionFitter(n_draws = 30L)
    gen <- forward_sim_generator(
      fitter = fitter,
      fit_spec = list(formula = y ~ x),
      pilot_bundle = pilot_bundle,
      predictor_generator = gaussian_predictors
    )
    bundle <- gen(list(n = 20L), list(task_id = "t1", rep_idx = 1L))

    expect_setequal(names(bundle$true_params), c("Intercept", "x", "sigma"))
    expect_setequal(bundle$vars_of_interest, c("Intercept", "x", "sigma"))
    expect_equal(names(bundle$true_params), bundle$vars_of_interest)
  })

  it("pins theta deterministically by rep_idx", {
    withr::local_seed(13)
    fitter <- LinearRegressionFitter(n_draws = 40L)
    gen <- forward_sim_generator(
      fitter = fitter,
      fit_spec = list(formula = y ~ x),
      pilot_bundle = pilot_bundle,
      predictor_generator = gaussian_predictors
    )
    b1 <- gen(list(n = 20L), list(task_id = "a", rep_idx = 1L))
    b2 <- gen(list(n = 20L), list(task_id = "b", rep_idx = 1L))
    expect_equal(b1$true_params, b2$true_params)

    b3 <- gen(list(n = 20L), list(task_id = "c", rep_idx = 2L))
    expect_false(identical(b1$true_params, b3$true_params))
    expect_equal(b3$meta$truth_draw_id, 2L)
  })

  it("wraps rep_idx modulo n_draws", {
    withr::local_seed(14)
    n_draws <- 25L
    fitter <- LinearRegressionFitter(n_draws = n_draws)
    gen <- forward_sim_generator(
      fitter = fitter,
      fit_spec = list(formula = y ~ x),
      pilot_bundle = pilot_bundle,
      predictor_generator = gaussian_predictors
    )
    rep_idx_wrap <- n_draws + 5L
    expected_id <- ((rep_idx_wrap - 1L) %% n_draws) + 1L
    bw <- gen(list(n = 20L), list(task_id = "tw", rep_idx = rep_idx_wrap))
    expect_equal(bw$meta$truth_draw_id, expected_id)
  })

  it("rejects a non-Fitter object", {
    expect_error(
      forward_sim_generator(
        fitter = "not_a_fitter",
        fit_spec = list(formula = y ~ x),
        pilot_bundle = pilot_bundle,
        predictor_generator = gaussian_predictors
      ),
      class = "bayesim_config_error"
    )
  })
})
