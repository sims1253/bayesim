library(bayesim)
library(extraDistr)
library(brms)
library(testthat)

# Unit tests for custom beta-prime

# Just a small function, to get the alpha shape-argument, for reference beta_prime functions
get_alpha <- function(mu, beta) {
  alpha <- mu * (beta - 1)
  return(alpha)
}

n <- 10000   # number of testvalues
eps <- 1e-6
x <- exp(seq(from = eps , to = 200 , length.out = n)) # testset, exp(200) comes close to Max-Double
unit <- seq(from = eps, to = 1 - eps, length.out = n)

n_small <- 10
mus <- seq(from = eps, to = 20, length.out = n_small)
betas <- seq(from = 1 + eps, to = 20, length.out = n_small)

test_that("custom-betaprime", {
  # calculate beta-prime
  dbetaprime_results <- bayesim::dbetaprime(x, mu = 1, beta = 2)
  qbetaprime_results <- bayesim::qbetaprime(unit, mu = 1, beta = 2)
  # check length
  expect_equal(n, length(dbetaprime_results))
  expect_equal(n, length(qbetaprime_results))
  # check values against comparable implementation
  expect_eps(dbetaprime_results, extraDistr::dbetapr(x, get_alpha(1, 2), 2), eps)
  expect_eps(qbetaprime_results, extraDistr::qbetapr(unit, get_alpha(1, 2), 2), eps)
  # also check other shape parameters
  expect_eps(bayesim::dbetaprime(x, mu = 1, beta = 4), extraDistr::dbetapr(x, get_alpha(1, 4), 4), eps)
  expect_eps(bayesim::qbetaprime(unit, mu = 1, beta = 4), extraDistr::qbetapr(unit, get_alpha(1, 4), 4), eps)

  # check many shape parameters
  for(m in mus) {
    for(b in betas) {
      expect_eps(bayesim::dbetaprime(x, mu = m, beta = b), extraDistr::dbetapr(x, get_alpha(m, b), b), eps)
      expect_eps(bayesim::qbetaprime(unit, mu = m, beta = b), extraDistr::qbetapr(unit, get_alpha(m, b), b), eps)
    }
  }

  # now check some errors
  expect_error(bayesim::dbetaprime(-1, mu = 2, beta = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dbetaprime(1, mu = 0, beta = 2))  # mu is not allowed to be 0 or smaller
  expect_error(bayesim::dbetaprime(1, mu = 1, beta = 1))  # beta is not allowed to be 1 or smaller
  expect_error(bayesim::qbetaprime(c(-1, 2), mu = 2, beta = 2))  # q is not allowed to be outside [0, 1]
  expect_error(bayesim::qbetaprime(1, mu = 0, beta = 2))  # mu is not allowed to be 0 or smaller
  expect_error(bayesim::qbetaprime(1, mu = 1, beta = 1))  # beta is not allowed to be 1 or smaller
})
