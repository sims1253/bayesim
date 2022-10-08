
test_that("test custom expect_eps function", {
  # wrong amount of arguments
  expect_error(expect_eps(1, 1))
  expect_error(expect_eps(1, 1, 1, 1, 1))
  # usage of non-numeric types, which is disallowed
  expect_error(expect_eps(c("a", 1), c("b", 1.1), c(2, 0.2)))
  # all scalars
  expect_success(expect_eps(1, 1.1, 0.2))
  expect_success(expect_eps(1, 3, 3))
  expect_success(expect_eps(-1, -1.1, 0.2))
  expect_failure(expect_eps(1, 2, 0.2))
  expect_failure(expect_eps(2, 1, 0.2))
  expect_error(expect_eps(0.1, 0.11, -0.1))
  # a vector, b scalar, eps scalar
  expect_success(expect_eps(c(1, 1), 1.1, 0.2))
  expect_success(expect_eps(c(2, 1, 1), 1.1, eps = 0.2, r = 0.4))
  expect_failure(expect_eps(c(1, 2), 1.1, 0.2))
  expect_failure(expect_eps(c(1, 1.1), 2, 0.2))
  expect_failure(expect_eps(c(2, 1, 1), 1.1, eps = 0.2, r = 0.2))
  # a vector, b vector, eps scalar
  expect_success(expect_eps(c(1, 1), c(1.1, 1.1), 0.2))
  expect_success(expect_eps(c(2, 1, 1), c(1.1, 1.1, 1.1), eps = 0.2, r = 0.4))
  expect_failure(expect_eps(c(1, 2), c(1.1, 1.1), 0.2))
  expect_failure(expect_eps(c(2, 1, 1), c(1.1, 1.1, 1.1), eps = 0.2, r = 0.2))
  expect_error(expect_eps(c(1, 1, 1), c(1.1, 1.1), 0.2))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1), 0.2))
  # a vector, b scalar, eps vector
  expect_success(expect_eps(c(1, 1), 1.1, c(0.2, 0.3)))
  expect_success(expect_eps(c(1, 1, 2), 1.1, eps = c(0.2, 0.3, 0.2), r = 0.4))
  expect_success(expect_eps(c(1, 1), 1.1, c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 2), 1.1, c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 1), 2, c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 1, 2), 1.1, eps = c(0.2, 0.3, 0.2), r = 0.2))
  expect_error(expect_eps(c(1, 1, 1), 1.1, c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 2, 3), 1.1, c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 1), 1.1, c(0.2, -0.3)))
  # all vectors
  expect_success(expect_eps(c(1, 1), c(1.1, 1.1), c(0.2, 0.3)))
  expect_success(expect_eps(c(1, 3), c(1.1, 3.1), c(0.2, 0.3)))
  expect_success(expect_eps(c(1, 5), c(1.1, 7), c(0.2, 3)))
  expect_success(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps = c(0.2, 0.3, 0.2), r = 0.4))
  expect_failure(expect_eps(c(1, 2), c(1.1, 1.1), c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps = c(0.2, 0.3, 0.2), r = 0.2))
  expect_error(expect_eps(c(1, 1, 1), c(1.1, 1.1), c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1), c(0.2, 0.3)))
  # additional r-error-tests
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps = c(0.2, 0.3, 0.2), r = -0.1))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps = c(0.2, 0.3, 0.2), r = 1.1))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps = c(0.2, 0.3, 0.2), r = c(0.1, 0.2, 0.3)))
  expect_error(expect_eps(1, 1.1, 0.2, "r"))
})

test_that("test the test_rng-wrapper", {
  # Test-RNG constants
  eps <- 1e-6
  n_small <- 10
  n <- 10000
  accepted_means_eps <- 0.35
  p_acceptable_failures <- 0.05
  mus <- seq(from = 1 + eps, to = 10, length.out = n_small)
  alphas_r <- seq(from = 2 + eps, to = 10, length.out = n_small)

  # the one test, if it works (not much point, in checking any other expected state)
  expect_success(test_rng(
    rng_fun = bayesim::rlomax, metric_mu = mean, n = n, mus = mus, shapes = alphas_r,
    mu_eps = accepted_means_eps, p_acceptable_failures = p_acceptable_failures
  ))
  # if the margins are too low, the function has a high likelyhood to fail (for eps=0, will always fail)
  expect_failure(test_rng(
    rng_fun = bayesim::rlomax, metric_mu = mean, n = n, mus = mus, shapes = alphas_r,
    mu_eps = 0, p_acceptable_failures = 0
  ))

  # else check all forbidden arguments (jay!)
  # non-function type function arguments
  expect_error(test_rng(0,
    metric_mu = mean, n = n, mus = mus, shapes = alphas_r,
    mu_eps = 0, p_acceptable_failures = 0
  ))
  expect_error(test_rng(
    rng_fun = bayesim::rlomax, 0, n = n, mus = mus, shapes = alphas_r,
    mu_eps = 0, p_acceptable_failures = 0
  ))
  # switched function arguments
  expect_error(test_rng(
    rng_fun = mean, metric_mu = bayesim::rlomax, n = n, mus = mus, shapes = alphas_r,
    mu_eps = 0, p_acceptable_failures = 0
  ))
  # non numeric sample number argument
  expect_error(test_rng(
    rng_fun = bayesim::rlomax, metric_mu = mean, n = "R", mus = mus, shapes = alphas_r,
    mu_eps = accepted_means_eps, p_acceptable_failures = p_acceptable_failures
  ))
  # vector of numbers to sample
  expect_error(test_rng(
    rng_fun = bayesim::rlomax, metric_mu = mean, n = c(42, 73), mus = mus, shapes = alphas_r,
    mu_eps = accepted_means_eps, p_acceptable_failures = p_acceptable_failures
  ))
  # non numeric eps argument
  expect_error(test_rng(
    rng_fun = bayesim::rlomax, metric_mu = mean, n = n, mus = mus, shapes = alphas_r,
    mu_eps = "R", p_acceptable_failures = p_acceptable_failures
  ))
  # negative eps argument
  expect_error(test_rng(
    rng_fun = bayesim::rlomax, metric_mu = mean, n = n, mus = mus, shapes = alphas_r,
    mu_eps = -1, p_acceptable_failures = p_acceptable_failures
  ))
  # vector eps argument
  expect_error(test_rng(
    rng_fun = bayesim::rlomax, metric_mu = mean, n = n, mus = mus, shapes = alphas_r,
    mu_eps = c(0.1, 0.2), p_acceptable_failures = p_acceptable_failures
  ))
  # non numeric p argument
  expect_error(test_rng(
    rng_fun = bayesim::rlomax, metric_mu = mean, n = n, mus = mus, shapes = alphas_r,
    mu_eps = accepted_means_eps, p_acceptable_failures = "R"
  ))
  # too small p argument (smaller 0)
  expect_error(test_rng(
    rng_fun = bayesim::rlomax, metric_mu = mean, n = n, mus = mus, shapes = alphas_r,
    mu_eps = accepted_means_eps, p_acceptable_failures = -1
  ))
  # too big p argument (bigger 1)
  expect_error(test_rng(
    rng_fun = bayesim::rlomax, metric_mu = mean, n = n, mus = mus, shapes = alphas_r,
    mu_eps = accepted_means_eps, p_acceptable_failures = 2
  ))
  # vector p argument
  expect_error(test_rng(
    rng_fun = bayesim::rlomax, metric_mu = mean, n = n, mus = mus, shapes = alphas_r,
    mu_eps = accepted_means_eps, p_acceptable_failures = c(0.1, 0.2)
  ))
})

test_that("test custom expect_bigger", {
  # tests, that should be correct
  expect_success(expect_bigger(1, 0))
  expect_success(expect_bigger(c(1, 2), 0))
  expect_success(expect_bigger(1, c(-1, 0)))
  expect_success(expect_bigger(c(1, 2), c(-1, 0)))
  # tests, that should be incorrect
  expect_failure(expect_bigger(0, 1))
  expect_failure(expect_bigger(0, c(1, 2)))
  expect_failure(expect_bigger(c(-1, 0), 1))
  expect_failure(expect_bigger(c(-1, 0), c(1, 2)))
  expect_failure(expect_bigger(1, c(0, 2)))
  expect_failure(expect_bigger(c(2, 0), c(1, 2)))
  # test, that are expected to throw an error
  expect_error(expect_bigger(1))
  expect_error(expect_bigger(1, 2, 3))
  expect_error(expect_bigger(c(1, 2, 3), c(-1, 0)))
  expect_error(expect_bigger("R", 1))
  expect_error(expect_bigger(c(), 1))
})

test_that("test_brms_quantile", {
  # To test this, first construct a BRMS model with Cloglognormal
  n_brms <- 1000
  intercept <- 0.3
  sigma <- 0.6
  tresh <- 0.05

  # save old seed, to reset it later
  old_seed <- .Random.seed
  # Set predefined seed. Generating correct and "random" RNG data is not part of the BRMS recovery test.
  set.seed(9001)
  cloglog_data <- bayesim::rcloglognormal(n_brms, intercept, sigma)
  set.seed(old_seed)
  # Now that the data was generated, reset the old seed (as if nothing ever happened)

  # limit the interval. Cloglognormal BRMS is very sensitive for data at the boundary.
  eps_brms <- 1e-12
  allowed_interval <- c(eps_brms, 1 - eps_brms)
  cloglog_data <- bayesim::limit_data(cloglog_data, allowed_interval)

  # special BRMS test implementation (as it uses a simplified y ~ 1 model)
  BBmisc::suppressAll({
    fit <- brms::brm(y ~ 1,
                     family = bayesim::cloglognormal(), stanvars = bayesim::cloglognormal()$stanvars,
                     backend = "cmdstanr", cores = 2, data = list(y = cloglog_data)
    )
  })

  # OK, after all that preamble, now it is getting interesting!
  expect_true(bayesim::test_brms_quantile(fit, "b_Intercept", intercept, tresh) &&
                bayesim::test_brms_quantile(fit, "sigma", sigma, tresh))
  # This test should be correct (is the almost the same as in the Cloglog Testthat)
  expect_warning(expect_false(bayesim::test_brms_quantile(fit, "alpha", sigma, tresh)))
  # No alpha in Cloglognormal, which return false and throws a warning
  sigma_data <- posterior::extract_variable_matrix(fit, variable = "sigma")
  median_sigma <- median(sigma_data)
  expect_false(bayesim::test_brms_quantile(fit, "sigma", 2*tresh + 2*median_sigma, tresh))
  # definitively data not within quantiles
  expect_error(bayesim::test_brms_quantile())
  # wrong amount of arguments
  expect_error(bayesim::test_brms_quantile(c(1 ,2, 3), "sigma", sigma, tresh))
  # vector is not type BRMS (which is a R list)
  expect_error(bayesim::test_brms_quantile(fit, sigma, sigma, tresh))
  # sigma is not a string
  expect_error(bayesim::test_brms_quantile(fit, c("sigma", "b_Intercept"), sigma, tresh))
  # only one value at a time
  expect_error(bayesim::test_brms_quantile(fit, "sigma", "sigma", tresh))
  # reference value should be real scalar, not string
  expect_error(bayesim::test_brms_quantile(fit, "sigma", NA, tresh))
  # reference value may not be NA
  expect_error(bayesim::test_brms_quantile(fit, "sigma", c(sigma, b_Intercept), tresh))
  # reference value should not be a vector
  expect_error(bayesim::test_brms_quantile(fit, "sigma", sigma, "tresh"))
  # threshold has to be of type real
  expect_error(bayesim::test_brms_quantile(fit, "sigma", sigma, -1))
  # threshold has to be in the unit-interval
  expect_error(bayesim::test_brms_quantile(fit, "sigma", sigma, c(0.1, 0.5, 0.9)))
  # threshold has to have 2 entries at most
  expect_error(bayesim::test_brms_quantile(fit, "sigma", sigma, c()))
  # no entries for threshold is forbidden as well
  expect_error(bayesim::test_brms_quantile(fit, "sigma", sigma, NA))
  # threshold cannot be NA
  expect_error(bayesim::test_brms_quantile(fit, "sigma", sigma, 0.6))
  # for 0.6, the threshold would create bounds, with lowerbound > upperbound
  expect_error(bayesim::test_brms_quantile(fit, "sigma", sigma, thresh, debug = "TRUE"))
  # debug has to be of type boolean
  expect_error(bayesim::test_brms_quantile(fit, "sigma", sigma, thresh, debug = c(TRUE, FALSE)))
  # debug has to be a single boolean
  expect_error(bayesim::test_brms_quantile(fit, "sigma", sigma, thresh, debug = 0))
  # debug has to be type boolean

})
