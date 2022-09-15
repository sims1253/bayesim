
test_that("test custom expect_eps function", {
  # wrong amount of arguments
  expect_error(expect_eps(1, 1))
  expect_error(expect_eps(1, 1, 1, 1, 1))
  # usage of non-numeric types, which is disallowed
  expect_error(expect_eps(c("a", 1), c("b", 1.1), c(2, 0.2)))
  # all scalars
  expect_success(expect_eps(1, 1.1, 0.2))
  expect_success(expect_eps(1, 3, 3))
  expect_success(expect_eps(-1, -1.1, 0.2))
  expect_failure(expect_eps(1, 2, 0.2))
  expect_failure(expect_eps(2, 1, 0.2))
  expect_error(expect_eps(0.1, 0.11, -0.1))
  # a vector, b scalar, eps scalar
  expect_success(expect_eps(c(1, 1), 1.1, 0.2))
  expect_success(expect_eps(c(2, 1, 1), 1.1, eps = 0.2, r = 0.4))
  expect_failure(expect_eps(c(1, 2), 1.1, 0.2))
  expect_failure(expect_eps(c(1, 1.1), 2, 0.2))
  expect_failure(expect_eps(c(2, 1, 1), 1.1, eps = 0.2, r = 0.2))
  # a vector, b vector, eps scalar
  expect_success(expect_eps(c(1, 1), c(1.1, 1.1), 0.2))
  expect_success(expect_eps(c(2, 1, 1), c(1.1, 1.1, 1.1), eps = 0.2, r = 0.4))
  expect_failure(expect_eps(c(1, 2), c(1.1, 1.1), 0.2))
  expect_failure(expect_eps(c(2, 1, 1), c(1.1, 1.1, 1.1), eps = 0.2, r = 0.2))
  expect_error(expect_eps(c(1, 1, 1), c(1.1, 1.1), 0.2))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1), 0.2))
  # a vector, b scalar, eps vector
  expect_success(expect_eps(c(1, 1), 1.1, c(0.2, 0.3)))
  expect_success(expect_eps(c(1, 1, 2), 1.1, eps = c(0.2, 0.3, 0.2), r = 0.4))
  expect_success(expect_eps(c(1, 1), 1.1, c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 2), 1.1, c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 1), 2, c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 1, 2), 1.1, eps = c(0.2, 0.3, 0.2), r = 0.2))
  expect_error(expect_eps(c(1, 1, 1), 1.1, c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 2, 3), 1.1, c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 1), 1.1, c(0.2, -0.3)))
  # all vectors
  expect_success(expect_eps(c(1, 1), c(1.1, 1.1), c(0.2, 0.3)))
  expect_success(expect_eps(c(1, 3), c(1.1, 3.1), c(0.2, 0.3)))
  expect_success(expect_eps(c(1, 5), c(1.1, 7), c(0.2, 3)))
  expect_success(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps = c(0.2, 0.3, 0.2), r = 0.4))
  expect_failure(expect_eps(c(1, 2), c(1.1, 1.1), c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps = c(0.2, 0.3, 0.2), r = 0.2))
  expect_error(expect_eps(c(1, 1, 1), c(1.1, 1.1), c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1), c(0.2, 0.3)))
  # additional r-error-tests
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps = c(0.2, 0.3, 0.2), r = -0.1))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps = c(0.2, 0.3, 0.2), r = 1.1))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps = c(0.2, 0.3, 0.2), r = c(0.1, 0.2, 0.3)))
  expect_error(expect_eps(1, 1.1, 0.2, "r"))
})

test_that("test the test_rng-wrapper", {
  # Test-RNG constants
  eps <- 1e-6
  n_small <- 10
  n <- 10000
  accepted_means_eps <- 0.3
  p_acceptable_failures <- 0.05
  mus <- seq(from = 1 + eps, to = 10, length.out = n_small)
  alphas_r <- seq(from = 2 + eps, to = 10, length.out = n_small)

  # the one test, if it works (not much point, in checking any other expected state)
  expect_success(test_rng(rng_fun=bayesim::rlomax, metric_mu=mean, n=n, mus=mus, shapes=alphas_r,
           mu_eps=accepted_means_eps, p_acceptable_failures=p_acceptable_failures))
  # if the margins are too low, the function has a high likelyhood to fail (for eps=0, will always fail)
  expect_failure(test_rng(rng_fun=bayesim::rlomax, metric_mu=mean, n=n, mus=mus, shapes=alphas_r,
                          mu_eps=0, p_acceptable_failures=0))

  # else check all forbidden arguments (jay!)
  # non-function type function arguments
  expect_error(test_rng(0, metric_mu=mean, n=n, mus=mus, shapes=alphas_r,
                        mu_eps=0, p_acceptable_failures=0))
  expect_error(test_rng(rng_fun=bayesim::rlomax, 0, n=n, mus=mus, shapes=alphas_r,
                        mu_eps=0, p_acceptable_failures=0))
  # switched function arguments
  expect_error(test_rng(rng_fun=mean, metric_mu=bayesim::rlomax, n=n, mus=mus, shapes=alphas_r,
                        mu_eps=0, p_acceptable_failures=0))
  # non numeric sample number argument
  expect_error(test_rng(rng_fun=bayesim::rlomax, metric_mu=mean, n="R", mus=mus, shapes=alphas_r,
                          mu_eps=accepted_means_eps, p_acceptable_failures=p_acceptable_failures))
  # vector of numbers to sample
  expect_error(test_rng(rng_fun=bayesim::rlomax, metric_mu=mean, n=c(42, 73), mus=mus, shapes=alphas_r,
                        mu_eps=accepted_means_eps, p_acceptable_failures=p_acceptable_failures))
  # non numeric eps argument
  expect_error(test_rng(rng_fun=bayesim::rlomax, metric_mu=mean, n=n, mus=mus, shapes=alphas_r,
                        mu_eps="R", p_acceptable_failures=p_acceptable_failures))
  # negative eps argument
  expect_error(test_rng(rng_fun=bayesim::rlomax, metric_mu=mean, n=n, mus=mus, shapes=alphas_r,
                        mu_eps=-1, p_acceptable_failures=p_acceptable_failures))
  # vector eps argument
  expect_error(test_rng(rng_fun=bayesim::rlomax, metric_mu=mean, n=n, mus=mus, shapes=alphas_r,
                        mu_eps=c(0.1, 0.2), p_acceptable_failures=p_acceptable_failures))
  # non numeric p argument
  expect_error(test_rng(rng_fun=bayesim::rlomax, metric_mu=mean, n=n, mus=mus, shapes=alphas_r,
                        mu_eps=accepted_means_eps, p_acceptable_failures="R"))
  # too small p argument (smaller 0)
  expect_error(test_rng(rng_fun=bayesim::rlomax, metric_mu=mean, n=n, mus=mus, shapes=alphas_r,
                        mu_eps=accepted_means_eps, p_acceptable_failures=-1))
  # too big p argument (bigger 1)
  expect_error(test_rng(rng_fun=bayesim::rlomax, metric_mu=mean, n=n, mus=mus, shapes=alphas_r,
                        mu_eps=accepted_means_eps, p_acceptable_failures=2))
  # vector p argument
  expect_error(test_rng(rng_fun=bayesim::rlomax, metric_mu=mean, n=n, mus=mus, shapes=alphas_r,
                        mu_eps=accepted_means_eps, p_acceptable_failures=c(0.1, 0.2)))
})

test_that("test custom expect_bigger", {
  # tests, that should be correct
  expect_success(expect_bigger(1, 0))
  expect_success(expect_bigger(c(1, 2), 0))
  expect_success(expect_bigger(1, c(-1, 0)))
  expect_success(expect_bigger(c(1, 2), c(-1, 0)))
  # tests, that should be incorrect
  expect_failure(expect_bigger(0, 1))
  expect_failure(expect_bigger(0, c(1, 2)))
  expect_failure(expect_bigger(c(-1, 0), 1))
  expect_failure(expect_bigger(c(-1, 0), c(1, 2)))
  expect_failure(expect_bigger(1, c(0, 2)))
  expect_failure(expect_bigger(c(2, 0), c(1, 2)))
  # test, that are expected to throw an error
  expect_error(expect_bigger(1))
  expect_error(expect_bigger(1, 2, 3))
  expect_error(expect_bigger(c(1, 2, 3), c(-1, 0)))
  expect_error(expect_bigger("R", 1))
  expect_error(expect_bigger(c(), 1))
})

test_that("test and document the rest of test-helper.R", {
  skip("test and document the rest of test-helper.R")
})
