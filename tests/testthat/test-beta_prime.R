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

n <- 100   # number of testvalues
x <- seq(from = 0 , to = 10 , length.out=n) # testset
eps <- 1e-3
test_that("dbetaprime", {
  expect_eps(bayesim::dbetaprime(x, mu = 1, beta = 2), extraDistr::dbetapr(x, get_alpha(1, 2), 2), eps)
})
