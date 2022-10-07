
test_that("custom-cloglognormal", {

  # load in values
  data <- readRDS("precalc_values/cloglognormal_refdata")
  pdf_data <- readRDS("precalc_values/cloglognormal_refpdf")

  # read data into variables. Makes it more readible and prevents rewriting code
  n <- data$n
  n_small <- data$n_small
  eps <- data$eps
  mus <- data$mus
  sigmas <- data$shapes

  x <- data$x
  pdf_ref <- as.matrix(pdf_data)

  # calculate beta-prime
  dcloglognormal_results <- bayesim::dcloglognormal(x, mu = 1, sigma = 2)
  # check length
  expect_equal(n, length(dcloglognormal_results))
  # check against one precalculated value
  expect_eps(0.4557343, bayesim::dcloglognormal(x = 0.5, mu = 1, sigma = 2), eps)

  # check the RNG will return the correct number of samples
  cloglognormal_samples <- bayesim::rcloglognormal(n, 2, 3)
  expect_equal(n, length(cloglognormal_samples))

  # check PDF data against precalculated reference data
  for (outer in 1:n_small) {
    for (inner in 1:n_small) {
      mu <- mus[outer]
      sigma <- sigmas[inner]
      expect_eps(bayesim::dcloglognormal(x, mu, sigma), pdf_ref[[outer, inner]], eps)
    }
  }

  n_rng <- 100000
  accepted_medians_eps <- 0.12
  p_acceptable_failures <- 0.05

  # check the RNG is not too far of the input value
  test_rng(
    rng_fun = bayesim::rcloglognormal, metric_mu = median, n = n_rng, mus = mus, shapes = sigmas,
    mu_eps = accepted_medians_eps, p_acceptable_failures = p_acceptable_failures, mu_link = cloglog
  )



  # now check density function for some errors
  expect_error(bayesim::dcloglognormal(0.5, 2)) # to few arguments
  expect_error(bayesim::dcloglognormal(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dcloglognormal(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dcloglognormal(2, mu = 2, sigma = 2)) # x is not allowed to be bigger 1
  expect_error(bayesim::dcloglognormal(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  expect_error(bayesim::dcloglognormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed


  # do same for RNG function
  expect_error(bayesim::rcloglognormal(100, 2)) # to few arguments
  expect_error(bayesim::rcloglognormal(10, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rcloglognormal(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rcloglognormal("r", mu = 2, sigma = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(bayesim::rcloglognormal(100, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller

  n_brms <- 1000
  intercept <- 0.5
  sigma <- 1
  thresh <- 0.05

  # save old seed, to reset it later
  old_seed <- .Random.seed
  # Set predefined seed. Generating correct and "random" RNG data is not part of the BRMS recovery test.
  set.seed(9001)
  cloglog_data <- bayesim::rcloglognormal(n_brms, intercept, sigma)
  cloglog_data_unchanged <- cloglog_data
  set.seed(old_seed)
  # Now that the data was generated, reset the old seed (as if nothing ever happened)

  # limit the interval. Cloglognormal BRMS is very sensitive for data at the boundary.
  eps_brms <- 1e-12
  allowed_interval <- c(eps_brms, 1 - eps_brms)
  cloglog_data_limited <- bayesim::limit_data(cloglog_data, allowed_interval)
  interval_str <- paste0("[", eps_brms, ", 1 - (", eps_brms, ")]")
  warning(paste0("Cloglog BRMS test with only simple model y ~ 1."))

  # BRMS clolognormal with limited data
  BBmisc::suppressAll({
    fit1 <- brms::brm(y ~ 1,
      family = bayesim::cloglognormal(), stanvars = bayesim::cloglognormal()$stanvars,
      backend = "cmdstanr", cores = 2, data = list(y = cloglog_data_limited)
    )
  })
  expect_true(bayesim::test_brms_quantile(fit1, "b_Intercept", intercept, thresh) &&
    bayesim::test_brms_quantile(fit1, "sigma", sigma, thresh))

  # BRMS clolognormal_limited (which will limit data itself) with not limited data
  BBmisc::suppressAll({
    fit2 <- brms::brm(y ~ 1,
                     family = bayesim::cloglognormal_limited(), stanvars = bayesim::cloglognormal_limited()$stanvars,
                     backend = "cmdstanr", cores = 2, data = list(y = cloglog_data_unchanged)
    )
  })
  expect_true(bayesim::test_brms_quantile(fit2, "b_Intercept", intercept, thresh) &&
                bayesim::test_brms_quantile(fit2, "sigma", sigma, thresh))

  # BRMS cloglognormal not with limited data, should fail (due to data too close to the boundaries)
  expect_error(
    BBmisc::suppressAll({
      fit3 <- brms::brm(y ~ 1,
                       family = bayesim::cloglognormal(), stanvars = bayesim::cloglognormal()$stanvars,
                       backend = "cmdstanr", cores = 2, data = list(y = cloglog_data_unchanged)
      )
    })
  )

  warning("Cloglognormal requires data to be limited (to around [1e-12, 1 - 1e-12]).
          Cloglognormal_limited should adress this issue (by limiting it in Stan code.
          Maximilian should check it sometime! :)")
})
