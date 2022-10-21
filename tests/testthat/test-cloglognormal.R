
test_that("custom-cloglognormal", {

  # load in values
  pdf_data <- readRDS("precalc_values/cloglognormal_refpdf")
  # intervals
  eps <- 1e-12 # 2 digits more, than sim
  unit_int <- c(eps, 1 - eps)
  mu_unit_int <- c(0.1, 0.9)
  pos_int <- c(eps, 200)
  shape_int <- c(0.1, 20)
  n <- 1000
  n_small <- 20

  mus <- cloglog(seq(from=mu_unit_int[1], to=mu_unit_int[2], length.out=n_small))
  sigmas <- seq(from=shape_int[1], to=shape_int[2], length.out=n_small)
  x <- seq(from=unit_int[1], to=unit_int[2], length.out=n)

  pdf_ref <- as.matrix(pdf_data)

  # calculate beta-prime
  dcloglognormal_results <- dcloglognormal(x, mu = 1, sigma = 2)
  # check length
  expect_equal(n, length(dcloglognormal_results))
  # check against one precalculated value
  expect_eps(0.4557343, dcloglognormal(x = 0.5, mu = 1, sigma = 2), eps=1e-6)

  # check the RNG will return the correct number of samples
  cloglognormal_samples <- rcloglognormal(n, 2, 3)
  expect_equal(n, length(cloglognormal_samples))

  # check PDF data against precalculated reference data
  for (outer in 1:n_small) {
    for (inner in 1:n_small) {
      mu <- mus[outer]
      sigma <- sigmas[inner]
      expect_eps(dcloglognormal(x, mu, sigma), pdf_ref[[outer, inner]], eps, relative=TRUE)
    }
  }

  n_rng <- 100000
  accepted_medians_eps <- 0.3
  p_acceptable_failures <- 0.1

  # check the RNG is not too far of the input value
  test_rng(
    rng_fun = rcloglognormal, metric_mu = median, n = n_rng, mu_list = mus, aux_par = sigmas,
    mu_eps = accepted_medians_eps, p_acceptable_failures = p_acceptable_failures, mu_link = cloglog
  )



  # now check density function for some errors
  expect_error(dcloglognormal(0.5, 2)) # to few arguments
  expect_error(dcloglognormal(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(dcloglognormal(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(dcloglognormal(2, mu = 2, sigma = 2)) # x is not allowed to be bigger 1
  expect_error(dcloglognormal(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  expect_error(dcloglognormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed


  # do same for RNG function
  expect_error(rcloglognormal(100, 2)) # to few arguments
  expect_error(rcloglognormal(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rcloglognormal(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(rcloglognormal("r", mu = 2, sigma = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(rcloglognormal(100, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller

  n_brms <- 1000
  intercept <- 0.5
  sigma <- 1
  thresh <- 0.05

  # save old seed, to reset it later
  old_seed <- .Random.seed
  # Set predefined seed. Generating correct and "random" RNG data is not part of the BRMS recovery test.
  set.seed(9001)
  cloglog_data <- rcloglognormal(n_brms, intercept, sigma)
  set.seed(old_seed)
  # Now that the data was generated, reset the old seed (as if nothing ever happened)

  # limit the interval. Cloglognormal BRMS is very sensitive for data at the boundary.
  eps_brms <- 1e-12
  allowed_interval <- c(eps_brms, 1 - eps_brms)
  cloglog_data <- limit_data(cloglog_data, allowed_interval)
  interval_str <- paste0("[", eps_brms, ", 1 - (", eps_brms, ")]")
  warning(paste0("Cloglog BRMS test with only simple model y ~ 1. And also manually limited data to: ", interval_str, "."))

  # special BRMS test implementation (as it uses a simplified y ~ 1 model)
  BBmisc::suppressAll({
    fit <- brms::brm(y ~ 1,
      family = cloglognormal(), stanvars = cloglognormal()$stanvars,
      backend = "cmdstanr", cores = 2, data = list(y = cloglog_data)
    )
  })
  expect_true(test_brms_quantile(fit, "b_Intercept", intercept, thresh) &&
    test_brms_quantile(fit, "sigma", sigma, thresh))
})
