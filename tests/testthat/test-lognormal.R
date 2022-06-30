library(bayesim)
library(brms)
library(testthat)

n <- 10000
eps <- 1e-6

n_small <- 10
x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
mus <- seq(from = -5, to = 20, length.out = n_small)
# weirdly enough, the log(median(RNG)) yields Infs for mu >= 4 ?
sigmas <- seq(from = 1 + eps, to = 20, length.out = n_small)

accepted_medians_eps <- 0.5
p_acceptable_failures <- 0.05

test_that("custom-lognormal", {
  # calculate beta-prime
  dlognormal_results <- bayesim::dlognormal(x, mu = 1, sigma = 2)
  # check length
  expect_equal(n, length(dlognormal_results))
  # check against one precalculated value
  expect_eps(0.278794, bayesim::dlognormal(x=0.5, mu=1, sigma=2), eps)

  warning("Think about, how to check dlognormal against a reference implementation, or precalculated values.")

  # check the RNG will return the correct number of samples
  lognormal_samples <- bayesim::rlognormal(n, 2, 3)
  expect_equal(n, length(lognormal_samples))

  # check the RNG is not too far of the input value
  test_rng(rng_fun=bayesim::rlognormal, metric_mu=median, n=n, mus=mus, shapes=sigmas,
           mu_eps=accepted_medians_eps, p_acceptable_failures=p_acceptable_failures, mu_link=log)



  # now check density function for some errors
  expect_error(bayesim::dlognormal(0.5, 2)) # to few arguments
  expect_error(bayesim::dlognormal(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dlognormal(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dlognormal(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  expect_error(bayesim::dlognormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed


  # do same for RNG function
  expect_error(bayesim::rlognormal(100, 2)) # to few arguments
  expect_error(bayesim::rlognormal(10, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rlognormal(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rlognormal("r", mu = 2, sigma = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(bayesim::rlognormal(100, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
})
