library(bayesim)
library(brms)
library(testthat)


test_that("custom-lognormal", {

  # load in values
  data <- readRDS("precalc_values/lognormal_refdata")
  pdf_data <- readRDS("precalc_values/lognormal_refpdf")

  n <- data$n
  n_small <- data$n_small
  eps <- data$eps
  mus <- data$mus
  sigmas <- data$shapes

  x <- data$x
  pdf_ref <- as.matrix(pdf_data)

  # calculate beta-prime
  dlognormal_custom_results <- bayesim::dlognormal_custom(x, mu = 1, sigma = 2)
  # check length
  expect_equal(n, length(dlognormal_custom_results))
  # check against one precalculated value
  expect_eps(0.278794, bayesim::dlognormal_custom(x=0.5, mu=1, sigma=2), eps)

  # check the RNG will return the correct number of samples
  lognormal_custom_samples <- bayesim::rlognormal_custom(n, 2, 3)
  expect_equal(n, length(lognormal_custom_samples))

  for(outer in 1:n_small) {
    for(inner in 1:n_small) {
      mu <- mus[outer]
      sigma <- sigmas[inner]
      expect_eps(bayesim::dlognormal_custom(x, mu, sigma), pdf_ref[[outer, inner]], eps)
    }
  }

  n_rng <- 100000
  accepted_medians_eps <- 0.15
  p_acceptable_failures <- 0.05
  # check the RNG is not too far of the input value
  test_rng(rng_fun=bayesim::rlognormal_custom, metric_mu=median, n=n_rng, mus=mus, shapes=sigmas,
           mu_eps=accepted_medians_eps, p_acceptable_failures=p_acceptable_failures, mu_link=log)



  # now check density function for some errors
  expect_error(bayesim::dlognormal_custom(0.5, 2)) # to few arguments
  expect_error(bayesim::dlognormal_custom(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dlognormal_custom(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dlognormal_custom(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  expect_error(bayesim::dlognormal_custom("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed


  # do same for RNG function
  expect_error(bayesim::rlognormal_custom(100, 2)) # to few arguments
  expect_error(bayesim::rlognormal_custom(10, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rlognormal_custom(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rlognormal_custom("r", mu = 2, sigma = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(bayesim::rlognormal_custom(100, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller

  expect_brms_family(link=identity, family=bayesim::lognormal_custom, rng=bayesim::rlognormal_custom, shape_name="sigma")
})
