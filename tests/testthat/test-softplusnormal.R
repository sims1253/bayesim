library(bayesim)
library(brms)
library(testthat)

n <- 10000
eps <- 1e-6

n_small <- 10
x <- seq(from = eps, to = 200, length.out = n) # testset, exp(200) comes close to Max-Double
mus <- seq(from = -5, to = 20, length.out = n_small)
# weirdly enough, the softplus(median(RNG)) yields Infs for mu >= 4 ?
sigmas <- seq(from = 1 + eps, to = 20, length.out = n_small)

accepted_medians_eps <- 0.5
p_acceptable_failures <- 0.05

test_that("custom-softplusnormal", {
  # calculate beta-prime
  dsoftplusnormal_results <- bayesim::dsoftplusnormal(x, mu = 1, sigma = 2)
  # check length
  expect_equal(n, length(dsoftplusnormal_results))
  # check against one precalculated value
  expect_eps(0.7800716, bayesim::dsoftplusnormal(x=0.5, mu=1, sigma=2), eps)

  warning("Think about, how to check dsoftplusnormal against a reference implementation, or precalculated values.")

  # now check density function for some errors
  expect_error(bayesim::dsoftplusnormal(0.5, 2)) # to few arguments
  expect_error(bayesim::dsoftplusnormal(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dsoftplusnormal(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dsoftplusnormal(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  expect_error(bayesim::dsoftplusnormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed

  skip("Softplusnormal-RNG throws NaNs constantly. Repair Softplusnormal and remove skip!")

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

  # check the RNG is not too far of the input value
  test_rng(rng_fun=bayesim::rsoftplusnormal, metric_mu=median, n=n, mus=mus, shapes=sigmas,
           mu_eps=accepted_medians_eps, p_acceptable_failures=p_acceptable_failures, mu_link=softplus)
})
