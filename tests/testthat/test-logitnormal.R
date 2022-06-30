library(bayesim)
library(brms)
library(testthat)

n <- 10000
eps <- 1e-6

n_small <- 10
x <- seq(from = eps, to = 1 - eps, length.out = n)
mus <- seq(from = -5, to = 20, length.out = n_small)
# weirdly enough, the logit(median(RNG)) yields Infs for mu >= 4 ?
sigmas <- seq(from = 1 + eps, to = 20, length.out = n_small)

accepted_medians_eps <- 0.5
p_acceptable_failures <- 0.05

test_that("custom-logitnormal", {
  # calculate beta-prime
  dlogitnormal_results <- bayesim::dlogitnormal(x, mu = 1, sigma = 2)
  # check length
  expect_equal(n, length(dlogitnormal_results))
  # check against one precalculated value
  expect_eps(0.7041307, bayesim::dlogitnormal(x=0.5, mu=1, sigma=2), eps)

  warning("Think about, how to check dlogitnormal against a reference implementation, or precalculated values.")

  # check the RNG will return the correct number of samples
  logitnormal_samples <- bayesim::rlogitnormal(n, 2, 3)
  expect_equal(n, length(logitnormal_samples))

  # check the RNG is not too far of the input value
  test_rng(rng_fun=bayesim::rlogitnormal, metric_mu=median, n=n, mus=mus, shapes=sigmas,
           mu_eps=accepted_medians_eps, p_acceptable_failures=p_acceptable_failures, mu_link=logit)



  # now check density function for some errors
  expect_error(bayesim::dlogitnormal(0.5, 2)) # to few arguments
  expect_error(bayesim::dlogitnormal(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dlogitnormal(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dlogitnormal(2, mu = 2, sigma = 2)) # x is not allowed to be bigger 1
  expect_error(bayesim::dlogitnormal(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  expect_error(bayesim::dlogitnormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed


  # do same for RNG function
  expect_error(bayesim::rlogitnormal(100, 2)) # to few arguments
  expect_error(bayesim::rlogitnormal(10, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rlogitnormal(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rlogitnormal("r", mu = 2, sigma = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(bayesim::rlogitnormal(100, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
})
