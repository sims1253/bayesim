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

n <- 100000 # number of testvalues
eps <- 1e-6
x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
unit <- seq(from = eps, to = 1 - eps, length.out = n)

n_small <- 10
mus <- seq(from = eps, to = 10, length.out = n_small)
alphas <- seq(from = 1 + eps, to = 10, length.out = n_small)

# mus_r <- seq(from = 1 + eps, to = 10, length.out = n_small)
alphas_r <- seq(from = 1.6 + eps, to = 10, length.out = n_small)
accepted_means_eps <- 0.2
p_acceptable_failures <- 0.05

test_that("custom-lomax", {
  # calculate lomax
  dlomax_results <- bayesim::dlomax(x, mu = 8, alpha = 2)
  qlomax_results <- bayesim::qlomax(unit, mu = 8, alpha = 2)
  # check length
  expect_equal(n, length(dlomax_results))
  expect_equal(n, length(qlomax_results))

  # check many shape parameters on pdf and qdf
  for (alpha in alphas) {
    for (m in mus) {
      expect_eps(bayesim::dlomax(x, mu = m, alpha = alpha), extraDistr::dlomax(x, get_lambda(m, alpha), alpha), eps)
      expect_eps(bayesim::qlomax(unit, mu = m, alpha = alpha), extraDistr::qlomax(unit, get_lambda(m, alpha), alpha), eps)
    }
  }


  # check the RNG will return the correct number of samples
  lomax_samples <- bayesim::rlomax(n, 1, 3)
  expect_equal(n, length(lomax_samples))

  # shape variable -> bound gets instable RNG, arbitrary bound instead with alpha_r
  test_rng(rng_fun=bayesim::rlomax, metric_mu=mean, n=n, mus=mus, shapes=alphas_r,
           mu_eps=accepted_means_eps, p_acceptable_failures=p_acceptable_failures)
  # check the RNG is not too far of the input value

  # now check density function for some errors
  expect_error(bayesim::dlomax(1, 2)) # to few arguments
  expect_error(bayesim::dlomax(1, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dlomax(-1, mu = 2, alpha = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dlomax("r", mu = 2, alpha = 2)) # non-numeric arguments are disallowed
  expect_error(bayesim::dlomax(1, mu = 0, alpha = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::dlomax(1, mu = 1, alpha = 0)) # alpha is not allowed to be 1 or smaller

  # do same for quantile function
  expect_error(bayesim::qlomax(1, 2)) # to few arguments
  expect_error(bayesim::qlomax(1, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::qlomax("r", mu = 2, alpha = 2)) # non-numeric arguments are disallowed
  expect_error(bayesim::qlomax(1, mu = 0, alpha = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::qlomax(1, mu = 1, alpha = 0)) # alpha is not allowed to be 0 or smaller
  expect_error(bayesim::qlomax(c(-1, 2), mu = 2, alpha = 2)) # q is not allowed to be outside [0, 1]

  # do same for RNG function
  expect_error(bayesim::rlomax(100, 2)) # to few arguments
  expect_error(bayesim::rlomax(10, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rlomax(-1, mu = 2, alpha = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rlomax("r", mu = 2, alpha = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(bayesim::rlomax(100, mu = 0, alpha = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::rlomax(100, mu = 1, alpha = 0)) # alpha is not allowed to be 0 or smaller

  expect_brms_family(link=exp, family=bayesim::lomax, rng=bayesim::rlomax, shape_name="alpha")
})
