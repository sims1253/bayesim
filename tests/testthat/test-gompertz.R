
# Unit tests for custom gompertz

# Just a small function, to get the alpha shape-argument, for reference gompertz functions
get_a <- function(mu, beta) {
  # a of extraDistr
  a <- -(beta * log(0.5)) / (exp(mu * beta) - 1)
  return(a)
}

n <- 100000 # number of testvalues
eps <- 1e-6
x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
unit <- seq(from = eps, to = 1 - eps, length.out = n)

n_small <- 10
mus <- seq(from = eps, to = 10, length.out = n_small)
betas <- seq(from = eps, to = 10, length.out = n_small)

# mus_r <- seq(from = 1 + eps, to = 10, length.out = n_small)
betas_r <- seq(from = 0.5 + eps, to = 10, length.out = n_small)
# parameter for RNG test (going to close to 0 makes it a bit unstable)
accepted_median_eps <- 0.05
p_acceptable_failures <- 0.05

test_that("custom-gompertz", {
  # calculate gompertz
  dgompertz_results <- bayesim::dgompertz(x, mu = 10, b = 1)
  qgompertz_results <- bayesim::qgompertz(unit, mu = 10, b = 1)
  # check length
  expect_equal(n, length(dgompertz_results))
  expect_equal(n, length(qgompertz_results))

  # check many shape parameters
  accepted_median_eps <- 0.01
  for (m in mus) {
    for (b in betas) {
      expect_eps(bayesim::dgompertz(x, mu = m, b = b), extraDistr::dgompertz(x, get_a(m, b), b), eps)
      expect_eps(bayesim::qgompertz(unit, mu = m, b = b), extraDistr::qgompertz(unit, get_a(m, b), b), eps)
    }
  }

  # check the RNG will return the correct number of samples
  gompertz_samples <- bayesim::rgompertz(n, 1, 3)
  expect_equal(n, length(gompertz_samples))

  test_rng(
    rng_fun = bayesim::rgompertz, metric_mu = median, n = n, mu_list = mus, aux_par = betas_r,
    mu_eps = accepted_median_eps, p_acceptable_failures = p_acceptable_failures
  )
  # check the RNG is not too far of the input value

  # now check density function for some errors
  expect_error(bayesim::dgompertz(1, 2)) # to few arguments
  expect_error(bayesim::dgompertz(1, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dgompertz(-1, mu = 2, beta = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dgompertz(1, mu = 0, beta = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::dgompertz(1, mu = 1, beta = 0)) # beta is not allowed to be 0 or smaller
  expect_error(bayesim::dgompertz("r", mu = 2, beta = 2)) # non-numeric argumen
  # check values against comparable implementationts are disallowed

  # do same for quantile function
  expect_error(bayesim::qgompertz(1, 2)) # to few arguments
  expect_error(bayesim::qgompertz(1, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::qgompertz(1, mu = 0, beta = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::qgompertz(1, mu = 1, beta = 0)) # beta is not allowed to be 0 or smaller
  expect_error(bayesim::qgompertz(c(-1, 2), mu = 2, beta = 2)) # q is not allowed to be outside [0, 1]
  expect_error(bayesim::qgompertz("r", mu = 2, beta = 2)) # non-numeric arguments are disallowed


  # do same for RNG function
  expect_error(bayesim::rgompertz(100, 2)) # to few arguments
  expect_error(bayesim::rgompertz(10, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rgompertz(-1, mu = 2, b = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rgompertz("r", mu = 2, b = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(bayesim::rgompertz(100, mu = 0, beta = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::rgompertz(100, mu = 1, beta = -1)) # beta is not allowed to be 0 or smaller

  expect_brms_family(link = exp, family = bayesim::gompertz, rng = bayesim::rgompertz, shape_name = "beta")
})
