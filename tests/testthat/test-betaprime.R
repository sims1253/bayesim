
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

# Global setup of testing space
n <- 10000
eps <- 1e-12
x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
unit <- seq(from = eps, to = 1 - eps, length.out = n)
n_small <- 20
mu_list <- seq(from = eps, to = 1000, length.out = n_small)
phi_list <- seq(from = eps, to = 50, length.out = n_small)
accepted_relative_error <- 1e-12
accepted_mean_error <- 0.05
p_acceptable_failures <- 0.05
# parameter for RNG test (going to close to 0 makes it a bit unstable)

test_that("custom-betaprime", {
  # calculate beta-prime
  dbetaprime_results <- dbetaprime(x, mu = 1, phi = 2)
  qbetaprime_results <- qbetaprime(unit, mu = 1, phi = 2)
  # check length
  expect_equal(n, length(dbetaprime_results))
  expect_equal(n, length(qbetaprime_results))

  # Compare density and quantile functions to extraDistr
  for (m in mu_list) {
    for (phi in phi_list) {
      bayesim:::expect_eps(dbetaprime(x, mu = m, phi = phi),
        extraDistr::dbetapr(x, get_alpha(m, phi), get_beta(phi)),
        eps = accepted_relative_error,
        relative = TRUE
      )
      bayesim:::expect_eps(qbetaprime(unit, mu = m, phi = phi),
        extraDistr::qbetapr(unit, get_alpha(m, phi), get_beta(phi)),
        eps = accepted_relative_error,
        relative = TRUE
      )
    }
  }

  # check the RNG will return the correct number of samples
  betaprime_samples <- bayesim::rbetaprime(n, 2, 3)
  expect_equal(n, length(betaprime_samples))

  # check if the RNG is close enough to the true mean in most cases
  bayesim:::test_rng(
    rng_fun = rbetaprime,
    metric_mu = mean,
    n = 10 * n,
    mu_list = mu_list,
    aux_par = phi_list,
    mu_eps = accepted_mean_error,
    p_acceptable_failures = p_acceptable_failures,
    relative = TRUE
  )

  # check if the RNG can recover the quantiles
  bayesim:::test_rng_quantiles(
    rng_fun = rbetaprime,
    quantile_fun = qbetaprime,
    n = 10 * n,
    mu_list = mu_list,
    aux_par = phi_list,
    eps = accepted_mean_error,
    quantiles = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99),
    p_acceptable_failures = 0.1,
    relative = TRUE
  )

  # now check density function for some errors
  expect_error(dbetaprime(1, 2)) # to few arguments
  expect_error(dbetaprime(1, 2, 3, 4, 5)) # to many arguments
  expect_error(dbetaprime(-1, mu = 2, phi = 2)) # x is not allowed to be smaller 0
  expect_error(dbetaprime(1, mu = 0, phi = 2)) # mu is not allowed to be 0 or smaller
  expect_error(dbetaprime(1, mu = 1, phi = 0)) # phi is not allowed to be 0 or smaller
  expect_error(dbetaprime("r", mu = 2, phi = 2)) # non-numeric arguments are disallowed

  # do same for quantile function
  expect_error(qbetaprime(1, 2)) # to few arguments
  expect_error(qbetaprime(1, 2, 3, 4, 5)) # to many arguments
  expect_error(qbetaprime(1, mu = 0, phi = 2)) # mu is not allowed to be 0 or smaller
  expect_error(qbetaprime(1, mu = 1, phi = 0)) # phi is not allowed to be 0 or smaller
  expect_error(qbetaprime(c(-1, 2), mu = 2, phi = 2)) # q is not allowed to be outside [0, 1]
  expect_error(qbetaprime("r", mu = 2, phi = 2)) # non-numeric arguments are disallowed


  # do same for RNG function
  expect_error(rbetaprime(100, 2)) # to few arguments
  expect_error(rbetaprime(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rbetaprime(-1, mu = 2, phi = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(rbetaprime("r", mu = 2, phi = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(rbetaprime(100, mu = 0, phi = 2)) # mu is not allowed to be 0 or smaller
  expect_error(rbetaprime(100, mu = 1, phi = 0)) # phi is not allowed to be 0 or smaller

  expect_brms_family(link = exp, family = betaprime, rng = rbetaprime, shape_name = "phi")
})
