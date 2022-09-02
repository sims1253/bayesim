
n <- 100000 # number of testvalues
eps <- 1e-3
x <- seq(from = eps, to = 1 - eps, length.out = n)

n_small <- 10
mus <- seq(from = eps, to = 1 - eps, length.out = n_small)
sigmas <- seq(from = eps, to = 5, length.out = n_small)

# Smaller n, given rsimplex is quite slow!
n_r_small <- 3
mus_r <- seq(from = eps, to = 1 - eps, length.out = n_r_small)
sigmas_r <- seq(from = 0.1 + eps, to = 5, length.out = n_r_small)
accepted_medians_eps <- 0.02
p_acceptable_failures <- 0.05

test_that("custom-simplex", {


  # calculate simplex
  dsimplex_results <- bayesim::dsimplex(x, mu = 0.8, sigma = 2)
  # check length
  expect_equal(n, length(dsimplex_results))

  # check the RNG will return the correct number of samples
  simplex_samples <- bayesim::rsimplex(n, 0.8, 3)
  expect_equal(n, length(simplex_samples))

  # shape variable -> bound gets instable RNG, arbitrary bound instead with p_r
  test_rng(rng_fun=bayesim::rsimplex, metric_mu=median, n=n, mus=mus_r, shapes=sigmas_r,
           mu_eps=accepted_medians_eps, p_acceptable_failures=p_acceptable_failures)
  # check the RNG is not too far of the input value

  # check many shape parameters on pdf
  for (s in sigmas) {
    for (m in mus) {
      expect_eps(bayesim::dsimplex(x, mu = m, sigma = s), rmutil::dsimplex(x, m, s^2), eps)
    }
  }

  # now check density function for some errors
  expect_error(bayesim::dsimplex(1, 0.8)) # to few arguments
  expect_error(bayesim::dsimplex(1, 0.8, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dsimplex(-1, mu = 0.8, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dsimplex("r", mu = 0.8, sigma = 2)) # non-numeric arguments are disallowed
  expect_error(bayesim::dsimplex(1, mu = 0, sigma = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::dsimplex(1, mu = 1, sigma = 2)) # mu is not allowed to be 1 or bigger
  expect_error(bayesim::dsimplex(1, mu = 0.8, sigma = 0)) # p is not allowed to be 1 or smaller

  # do same for RNG function
  expect_error(bayesim::rsimplex(100, 0.8)) # to few arguments
  expect_error(bayesim::rsimplex(100.82, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rsimplex(-1, mu = 0.8, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rsimplex("r", mu = 0.8, sigma = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(bayesim::rsimplex(100, mu = 0, sigma = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::rsimplex(100, mu = 1, sigma = 2)) # mu is not allowed to be 1 or bigger
  expect_error(bayesim::rsimplex(100, mu = 0.8, sigma = 0)) # p is not allowed to be 0 or smaller

  expect_error(bayesim::rkuramaswamy(100, mu = 0.8, sigma = 1)) # simplex has to be spelled correctly!!!
  # small inside joke, given, there is a 50% chance, I misspelled it again. :P

  expect_brms_family(link=brms::inv_logit_scaled, family=bayesim::simplex, rng=bayesim::rsimplex, shape_name="sigma")

})

