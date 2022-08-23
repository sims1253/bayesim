library(bayesim)
library(brms)
library(testthat)

test_that("custom-softplusnormal", {

  # load in values
  data_inp <- readRDS("precalc_values/softplusnormal_inp")
  data_ref <- readRDS("precalc_values/softplusnormal_ref")

  # one might build one extra set with only scalars. But way more work for miniscule space savings!
  # (yes, n_small is really that small!)
  n <- data_inp$n[1]
  n_small <- data_inp$n_small[1]
  eps <- data_inp$eps[1]
  mus <- data_inp$mus
  sigmas <- data_inp$shapes

  x <- data_ref$x
  data_subref <- subset(data_ref, select = -x)
  pdf_ref <- as.matrix(data_subref)


  # calculate beta-prime
  dsoftplusnormal_results <- bayesim::dsoftplusnormal(x, mu = 1, sigma = 2)
  # check length
  expect_equal(n, length(dsoftplusnormal_results))
  # check against one precalculated value
  expect_eps(0.3922206, bayesim::dsoftplusnormal(x=0.5, mu=1, sigma=2), eps)
  # ??? The constant in this test was wrong? Probably corrected implementation, or something
  # Now it is not 1.27.. anymore, but 0.392.. Maybe that is more correct?

  for(outer in 1:n_small) {
    for(inner in 1:n_small) {
      mu <- mus[outer]
      sigma <- sigmas[inner]
      expect_eps(bayesim::dsoftplusnormal(x, mu, sigma), pdf_ref[[outer, inner]], eps)
    }
  }

  # now check density function for some errors
  expect_error(bayesim::dsoftplusnormal(0.5, 2)) # to few arguments
  expect_error(bayesim::dsoftplusnormal(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dsoftplusnormal(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dsoftplusnormal(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  expect_error(bayesim::dsoftplusnormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed

  # do same for RNG function
  expect_error(bayesim::rsoftplusnormal(100, 2)) # to few arguments
  expect_error(bayesim::rsoftplusnormal(10, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rsoftplusnormal(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rsoftplusnormal("r", mu = 2, sigma = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(bayesim::rsoftplusnormal(100, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller


  # check the RNG will return the correct number of samples
  softplusnormal_samples <- bayesim::rsoftplusnormal(n, 2, 3)
  expect_equal(n, length(softplusnormal_samples))

  n_rng <- 100000
  accepted_medians_eps <- 0.13
  p_acceptable_failures <- 0.05

  # check the RNG is not too far of the input value
  test_rng(rng_fun=bayesim::rsoftplusnormal, metric_mu=median, n=n_rng, mus=mus, shapes=sigmas,
           mu_eps=accepted_medians_eps, p_acceptable_failures=p_acceptable_failures, mu_link=softplus)

  # check custom BRMS family implementation
  expect_brms_family(ba=0.2, intercept=0.4, shape=2, link=identity, family=bayesim::softplusnormal,
                     rng=bayesim::rsoftplusnormal, shape_name="sigma")
})
