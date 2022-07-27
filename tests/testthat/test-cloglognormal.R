library(bayesim)
library(brms)
library(testthat)

n <- 100000
eps <- 1e-6

n_small <- 10
x <- seq(from = eps, to = 1 - eps, length.out = n)
mus <- seq(from = -20, to = 3, length.out = n_small)
# weirdly enough, the cloglog(median(RNG)) yields Infs for mu >= 4 ?
sigmas <- seq(from = 1 + eps, to = 20, length.out = n_small)

accepted_medians_eps <- 0.12
p_acceptable_failures <- 0.05

test_that("custom-cloglognormal", {
  # calculate beta-prime
  dcloglognormal_results <- bayesim::dcloglognormal(x, mu = 1, sigma = 2)
  # check length
  expect_equal(n, length(dcloglognormal_results))
  # check against one precalculated value
  expect_eps(0.4557343, bayesim::dcloglognormal(x=0.5, mu=1, sigma=2), eps)

  warning("Think about, how to check dcloglognormal against a reference implementation, or precalculated values.")

  # check the RNG will return the correct number of samples
  cloglognormal_samples <- bayesim::rcloglognormal(n, 2, 3)
  expect_equal(n, length(cloglognormal_samples))

  # check the RNG is not too far of the input value
  test_rng(rng_fun=bayesim::rcloglognormal, metric_mu=median, n=n, mus=mus, shapes=sigmas,
           mu_eps=accepted_medians_eps, p_acceptable_failures=p_acceptable_failures, mu_link=cloglog)



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

  warning("BRMS missing")
})
