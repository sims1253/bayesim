library(bayesim)
library(rmutil)
library(brms)
library(testthat)

# Unit tests for custom inversegaussian_custom

# Just lambda small function, to get the shape shape-argument, for reference inversegaussian_custom functions
get_lambda <- function(mu, shape) {
  lambda <- mu * (shape - 1)
  # lambda of extraDistr is inverted, compared to lambda in my pdf
  return(1 / lambda)
}

n <- 100000 # number of testvalues
eps <- 1e-6
x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
unit <- seq(from = eps, to = 1 - eps, length.out = n)

n_small <- 10
mus <- seq(from = eps, to = 10, length.out = n_small)
shapes <- seq(from = eps, to = 10, length.out = n_small)

# mus_r <- seq(from = 1 + eps, to = 10, length.out = n_small)
shapes_r <- seq(from = 0.5 + eps, to = 10, length.out = n_small)
accepted_means_eps <- 0.1
p_acceptable_failures <- 0.05

test_that("custom-inversegaussian_custom", {
  # calculate inversegaussian_custom
  dinversegaussian_custom_results <- bayesim::dinversegaussian_custom(x, mu = 8, shape = 2)
  # check length
  expect_equal(n, length(dinversegaussian_custom_results))



  # check the RNG will return the correct number of samples
  inversegaussian_custom_samples <- bayesim::rinversegaussian_custom(n, 1, 3)
  expect_equal(n, length(inversegaussian_custom_samples))

  # shape variable -> bound gets instable RNG, arbitrary bound instead with shape_r
  test_rng(rng_fun=bayesim::rinversegaussian_custom, metric_mu=mean, n=n, mus=mus, shapes=shapes_r,
           mu_eps=accepted_means_eps, p_acceptable_failures=p_acceptable_failures)
  # check the RNG is not too far of the input value

  # now check density function for some errors
  expect_error(bayesim::dinversegaussian_custom(1, 2)) # to few arguments
  expect_error(bayesim::dinversegaussian_custom(1, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dinversegaussian_custom(-1, mu = 2, shape = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dinversegaussian_custom("r", mu = 2, shape = 2)) # non-numeric arguments are disallowed
  expect_error(bayesim::dinversegaussian_custom(1, mu = 0, shape = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::dinversegaussian_custom(1, mu = 1, shape = 0)) # shape is not allowed to be 1 or smaller

  # do same for quantile function
  expect_error(bayesim::qinversegaussian_custom(1, 2)) # to few arguments
  expect_error(bayesim::qinversegaussian_custom(1, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::qinversegaussian_custom("r", mu = 2, shape = 2)) # non-numeric arguments are disallowed
  expect_error(bayesim::qinversegaussian_custom(1, mu = 0, shape = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::qinversegaussian_custom(1, mu = 1, shape = 0)) # shape is not allowed to be 0 or smaller
  expect_error(bayesim::qinversegaussian_custom(c(-1, 2), mu = 2, shape = 2)) # q is not allowed to be outside [0, 1]

  # do same for RNG function
  expect_error(bayesim::rinversegaussian_custom(100, 2)) # to few arguments
  expect_error(bayesim::rinversegaussian_custom(10, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rinversegaussian_custom(-1, mu = 2, shape = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rinversegaussian_custom("r", mu = 2, shape = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(bayesim::rinversegaussian_custom(100, mu = 0, shape = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::rinversegaussian_custom(100, mu = 1, shape = 0)) # shape is not allowed to be 0 or smaller


  expect_brms_family(n=1000, ba=0.5, int=1, shape=2, link=exp,family=bayesim::inversegaussian_custom,
                     rng=bayesim::rinversegaussian_custom, shape_name="shape", thresh = 0.05)

  skip("Some issues (probably different parametrizations) between Bayesim and rmutil.")
  # check many shape parameters on pdf
  for (shape in shapes) {
    for (m in mus) {
      # check against BRMS implementation
      expect_eps(bayesim::dinversegaussian_custom(x, mu = m, shape = shape), rmutil::dinvgauss(x, m, shape), eps)
    }
  }

})
