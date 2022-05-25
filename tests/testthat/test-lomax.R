library(bayesim)
library(extraDistr)
library(brms)
library(testthat)

# Unit tests for custom lomax

# Just lambda small function, to get the alpha shape-argument, for reference lomax functions
get_lambda <- function(mu, alpha) {
  lambda <- mu * (alpha - 1)
  # lambda of extraDistr is inverted, compared to lambda in my pdf
  return(1 / lambda)
}

n <- 10000 # number of testvalues
eps <- 1e-6
x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
unit <- seq(from = eps, to = 1 - eps, length.out = n)

n_small <- 10
mus <- seq(from = eps, to = 10, length.out = n_small)
alphas <- seq(from = 1 + eps, to = 10, length.out = n_small)
mus_r <- seq(from = 1 + eps, to = 10, length.out = n_small)
alphas_r <- seq(from = 2 + eps, to = 10, length.out = n_small)

p_acceptable_failures <- 0.2 # with arbitrary median_eps of 0.1, about 8-20% of medians will be outside that range

test_that("custom-lomax", {
  # calculate lomax
  dlomax_results <- bayesim::dlomax(x, mu = 8, alpha = 2)
  qlomax_results <- bayesim::qlomax(unit, mu = 8, alpha = 2)
  # check length
  expect_equal(n, length(dlomax_results))
  expect_equal(n, length(qlomax_results))
  # check values against comparable implementation
  expect_eps(dlomax_results, extraDistr::dlomax(x, get_lambda(8, 2), 2), eps)
  expect_eps(qlomax_results, extraDistr::qlomax(unit, get_lambda(8, 2), 2), eps)
  # also check other shape parameters
  expect_eps(bayesim::dlomax(x, mu = 2, alpha = 2), extraDistr::dlomax(x, get_lambda(2, 2), 2), eps)
  expect_eps(bayesim::qlomax(unit, mu = 2, alpha = 2), extraDistr::qlomax(unit, get_lambda(2, 2), 2), eps)

  # check the RNG for one test with mu = 1 and alpha = 3
  mu <- 1
  accepted_mean_eps <- 0.5
  lomax_samples <- bayesim::rlomax(n, mu, 3)
  expect_equal(n, length(lomax_samples))
  expect_eps(mean(lomax_samples), mu, accepted_mean_eps) # this test should work most of the time, but might fail sometimes

  # check many shape parameters on pdf and qdf
  for (alpha in alphas) {
    for (m in mus) {
      expect_eps(bayesim::dlomax(x, mu = m, alpha = alpha), extraDistr::dlomax(x, get_lambda(m, alpha), alpha), eps)
      expect_eps(bayesim::qlomax(unit, mu = m, alpha = alpha), extraDistr::qlomax(unit, get_lambda(m, alpha), alpha), eps)
    }
  }

  # check many shapre parameters on RNG. Uses values further away, from definition edge
  for (alpha in alphas_r) {
    mu_r_length <- length(mus_r)
    lomax_rng_means <- vector(length = mu_r_length)
    for (i in 0:mu_r_length) {
      lomax_rng_means[i] <- mean(bayesim::rlomax(n, mu = mus_r[i], alpha = alpha))
    }
    expect_eps(lomax_rng_means, mus_r, accepted_mean_eps, p_acceptable_failures)
  }

  # now check density function for some errors
  expect_error(bayesim::dlomax(1, 2)) # to few arguments
  expect_error(bayesim::dlomax(1, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dlomax(-1, mu = 2, alpha = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dlomax("r", mu = 2, alpha = 2)) # non-numeric arguments are disallowed
  expect_error(bayesim::dlomax(1, mu = 0, alpha = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::dlomax(1, mu = 1, alpha = 0)) # alpha is not allowed to be 1 or smaller
  expect_error(bayesim::dlomax(1, mu = -1, alpha = 2)) # mu is not allowed to be smaller than 0
  expect_error(bayesim::dlomax(1, mu = 1, alpha = 0)) # alpha is not allowed to be 1 or smaller

  # do same for quantile function
  expect_error(bayesim::qlomax(1, 2)) # to few arguments
  expect_error(bayesim::qlomax(1, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::qlomax(-1, mu = 2, alpha = 2)) # q is not allowed to be smaller 0
  expect_error(bayesim::qlomax("r", mu = 2, alpha = 2)) # non-numeric arguments are disallowed
  expect_error(bayesim::qlomax(1, mu = 0, alpha = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::qlomax(1, mu = 1, alpha = 0)) # alpha is not allowed to be 0 or smaller
  expect_error(bayesim::qlomax(c(-1, 2), mu = 2, alpha = 2)) # q is not allowed to be outside [0, 1]
  expect_error(bayesim::qlomax(1, mu = -1, alpha = 2)) # mu is not allowed to be smaller than 0
  expect_error(bayesim::qlomax(1, mu = 1, alpha = 0)) # alpha is not allowed to be 1 or smaller

  # do same for RNG function
  expect_error(bayesim::rlomax(100, 2)) # to few arguments
  expect_error(bayesim::rlomax(10, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rlomax(-1, mu = 2, alpha = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rlomax("r", mu = 2, alpha = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(bayesim::rlomax(100, mu = 0, alpha = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::rlomax(100, mu = 1, alpha = 0)) # alpha is not allowed to be 0 or smaller
  expect_error(bayesim::rlomax(100, mu = -1, alpha = 2)) # mu is not allowed to be smaller than 0
  expect_error(bayesim::rlomax(100, mu = 1, alpha = 0)) # alpha is not allowed to be 1 or smaller
})
