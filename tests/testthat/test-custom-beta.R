
get_a <- function(mu, phi) {
  return(mu * phi)
}

get_b <- function(mu, phi) {
  return((1 - mu) * phi)
}

n <- 100000 # number of testvalues
eps <- 1e-3
x <- seq(from = eps, to = 1 - eps, length.out = n)

n_small <- 10
mus <- seq(from = eps, to = 1 - eps, length.out = n_small)
phis <- seq(from = eps, to = 5, length.out = n_small)

test_that("test-custom-beta", {
  # calculate beta_custom
  dbeta_custom_results <- bayesim::dbeta_custom(x, mu = 0.8, phi = 2)
  expect_equal(n, length(dbeta_custom_results))

  for (phi in phis) {
    for (mu in mus) {
      expect_eps(
        bayesim::dbeta_custom(x, mu = mu, phi = phi),
        dbeta(x, get_a(mu, phi), get_b(mu, phi)),
        eps
        )
    }
  }

  # now check density function for some errors
  expect_error(bayesim::dbeta_custom(1, 0.8)) # to few arguments
  expect_error(bayesim::dbeta_custom(1, 0.8, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dbeta_custom(-1, mu = 0.8, phi = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dbeta_custom("r", mu = 0.8, phi = 2)) # non-numeric arguments are disallowed
  expect_error(bayesim::dbeta_custom(1, mu = 0, phi = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::dbeta_custom(1, mu = 0.8, phi = 2)) # mu is not allowed to be 1 or bigger
  expect_error(bayesim::dbeta_custom(1, mu = 0.8, phi = 0)) # p is not allowed to be 1 or smaller
})
