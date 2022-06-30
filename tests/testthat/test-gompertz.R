library(bayesim)
library(extraDistr)
library(brms)
library(testthat)

# Unit tests for custom gompertz

# Just a small function, to get the alpha shape-argument, for reference gompertz functions
get_b <- function(mu, eta) {
  b <- (1 / mu) * log((-1 / eta) * log(1 / 2) + 1)
  return(b)
}
get_a <- function(mu, eta) {
  # a of extraDistr
  a <- get_b(mu, eta) * eta
  return(a)
}

n <- 10000 # number of testvalues
eps <- 1e-6
x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
unit <- seq(from = eps, to = 1 - eps, length.out = n)

n_small <- 10
mus <- seq(from = eps, to = 10, length.out = n_small)
etas <- seq(from = eps, to = 10, length.out = n_small)

# mus_r <- seq(from = 1 + eps, to = 10, length.out = n_small)
# etas_r <- seq(from = 1 + eps, to = 10, length.out = n_small)
accepted_median_eps <- 0.25
p_acceptable_failures <- 0.05

test_that("custom-gompertz", {
  # calculate gompertz
  dgompertz_results <- bayesim::dgompertz(x, mu = 1, eta = 0.1)
  qgompertz_results <- bayesim::qgompertz(unit, mu = 1, eta = 0.1)
  # check length
  expect_equal(n, length(dgompertz_results))
  expect_equal(n, length(qgompertz_results))

  # check many shape parameters
  for (m in mus) {
    for (eta in etas) {
      expect_eps(bayesim::dgompertz(x, mu = m, eta = eta), extraDistr::dgompertz(x, get_a(m, eta), get_b(m, eta)), eps)
      expect_eps(bayesim::qgompertz(unit, mu = m, eta = eta), extraDistr::qgompertz(unit, get_a(m, eta), get_b(m, eta)), eps)
    }
  }

  # check the RNG will return the correct number of samples
  gompertz_samples <- bayesim::rgompertz(n, 1, 3)
  expect_equal(n, length(gompertz_samples))

  test_rng(rng_fun=bayesim::rgompertz, metric_mu=median, n=n, mus=mus, shapes=etas,
           mu_eps=accepted_median_eps, p_acceptable_failures=p_acceptable_failures)
  # check the RNG is not too far of the input value


  # now check density function for some errors
  expect_error(bayesim::dgompertz(1, 2)) # to few arguments
  expect_error(bayesim::dgompertz(1, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dgompertz(-1, mu = 2, eta = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dgompertz(1, mu = 0, eta = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::dgompertz(1, mu = 1, eta = 0)) # eta is not allowed to be 0 or smaller
  expect_error(bayesim::dgompertz("r", mu = 2, eta = 2)) # non-numeric argumen
  # check values against comparable implementationts are disallowed

  # do same for quantile function
  expect_error(bayesim::qgompertz(1, 2)) # to few arguments
  expect_error(bayesim::qgompertz(1, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::qgompertz(1, mu = 0, eta = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::qgompertz(1, mu = 1, eta = 0)) # eta is not allowed to be 0 or smaller
  expect_error(bayesim::qgompertz(c(-1, 2), mu = 2, eta = 2)) # q is not allowed to be outside [0, 1]
  expect_error(bayesim::qgompertz("r", mu = 2, eta = 2)) # non-numeric arguments are disallowed


  # do same for RNG function
  expect_error(bayesim::rgompertz(100, 2)) # to few arguments
  expect_error(bayesim::rgompertz(10, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rgompertz(-1, mu = 2, eta = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rgompertz("r", mu = 2, eta = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(bayesim::rgompertz(100, mu = 0, eta = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::rgompertz(100, mu = 1, eta = -1)) # eta is not allowed to be 0 or smaller
})
