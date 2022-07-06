library(bayesim)
library(extraDistr)
library(brms)
library(testthat)

# Unit tests for custom beta-prime

# Just a small function, to get the alpha shape-argument, for reference beta_prime functions
get_alpha <- function(mu, phi) {
  alpha <- mu * (phi + 1)
  return(alpha)
}
get_beta <- function(phi) {
  beta <- phi + 2
  return(beta)
}

n <- 10000 # number of testvalues
eps <- 1e-6
x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
unit <- seq(from = eps, to = 1 - eps, length.out = n)

n_small <- 10
mus <- seq(from = eps, to = 20, length.out = n_small)
phis <- seq(from = 1 + eps, to = 20, length.out = n_small)

accepted_means_eps <- 0.2
p_acceptable_failures <- 0.05
# mus_r <- seq(from = 1 + eps, to = 20, length.out = n_small)
# phis_r <- seq(from = 2 + eps, to = 20, length.out = n_small)

test_that("custom-betaprime", {
  # skip("Beta-Prime tests need to be updated for new parametrisation!")
  # calculate beta-prime
  dbetaprime_results <- bayesim::dbetaprime(x, mu = 1, phi = 2)
  qbetaprime_results <- bayesim::qbetaprime(unit, mu = 1, phi = 2)
  # check length
  expect_equal(n, length(dbetaprime_results))
  expect_equal(n, length(qbetaprime_results))

  # check many shape parameters
  for (m in mus) {
    for (phi in phis) {
      expect_eps(bayesim::dbetaprime(x, mu = m, phi = phi), extraDistr::dbetapr(x, get_alpha(m, phi), get_beta(phi)), eps)
      expect_eps(bayesim::qbetaprime(unit, mu = m, phi = phi), extraDistr::qbetapr(unit, get_alpha(m, phi), get_beta(phi)), eps)
    }
  }

<<<<<<< HEAD:tests/testthat/test-betaprime.R
  # check the RNG will return the correct number of samples
  betaprime_samples <- bayesim::rbetaprime(n, 2, 3)
=======
  # check RNG
  n <- 100000
  mu <- 2
  betaprime_samples <- bayesim::rbetaprime(n, mu, 3)
>>>>>>> master:tests/testthat/test-beta_prime.R
  expect_equal(n, length(betaprime_samples))

  # check the RNG is not too far of the input value
  test_rng(rng_fun=bayesim::rbetaprime, metric_mu=mean, n=n, mus=mus, shapes=phis,
           mu_eps=accepted_means_eps, p_acceptable_failures=p_acceptable_failures)



  # now check density function for some errors
  expect_error(bayesim::dbetaprime(1, 2)) # to few arguments
  expect_error(bayesim::dbetaprime(1, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dbetaprime(-1, mu = 2, phi = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dbetaprime(1, mu = 0, phi = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::dbetaprime(1, mu = 1, phi = 0)) # phi is not allowed to be 0 or smaller
  expect_error(bayesim::dbetaprime("r", mu = 2, phi = 2)) # non-numeric arguments are disallowed

  # do same for quantile function
  expect_error(bayesim::qbetaprime(1, 2)) # to few arguments
  expect_error(bayesim::qbetaprime(1, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::qbetaprime(1, mu = 0, phi = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::qbetaprime(1, mu = 1, phi = 0)) # phi is not allowed to be 0 or smaller
  expect_error(bayesim::qbetaprime(c(-1, 2), mu = 2, phi = 2)) # q is not allowed to be outside [0, 1]
  expect_error(bayesim::qbetaprime("r", mu = 2, phi = 2)) # non-numeric arguments are disallowed


  # do same for RNG function
  expect_error(bayesim::rbetaprime(100, 2)) # to few arguments
  expect_error(bayesim::rbetaprime(10, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rbetaprime(-1, mu = 2, phi = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rbetaprime("r", mu = 2, phi = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(bayesim::rbetaprime(100, mu = 0, phi = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::rbetaprime(100, mu = 1, phi = 0)) # phi is not allowed to be 0 or smaller
})
