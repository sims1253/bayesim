library(bayesim)
library(VGAM)
library(testthat)

n <- 100000 # number of testvalues
eps <- 1e-6
x <- seq(from = eps, to = 1 - eps, length.out = n)

n_small <- 10
mus <- seq(from = eps, to = 1 - eps, length.out = n_small)
ps <- seq(from = eps, to = 10, length.out = n_small)

# mus_r <- seq(from = 1 + eps, to = 10, length.out = n_small)
ps_r <- seq(from = 0.6 + eps, to = 10, length.out = n_small)
accepted_means_eps <- 0.02
p_acceptable_failures <- 0.05

get_q <- function(mu, p) {
  value <- -log(2) / log1p(-mu^p)
  return(value)
}

test_that("custom-kumaraswamy", {
  # calculate kumaraswamy
  dkumaraswamy_results <- bayesim::dkumaraswamy(x, mu = 0.8, p = 2)
  qkumaraswamy_results <- bayesim::qkumaraswamy(x, mu = 0.8, p = 2)
  # check length
  expect_equal(n, length(dkumaraswamy_results))
  expect_equal(n, length(qkumaraswamy_results))

  # check many shape parameters on pdf and qdf
  for (p in ps) {
    for (m in mus) {
      expect_eps(bayesim::dkumaraswamy(x, mu = m, p = p), VGAM::dkumar(x, p, get_q(m, p)), eps)
      #expect_eps(bayesim::qkumaraswamy(x, mu = m, p = p), VGAM::qkumar(x, p, get_q(m, p)), eps)
      diff <- abs(bayesim::qkumaraswamy(x, mu = m, p = p) - VGAM::qkumar(x, p, get_q(m, p)))
      if(isTRUE(any(diff > eps))) {
        print(c("p, m", p, m))
        data <- VGAM::qkumar(x, p, get_q(m, p))
        print(data)
        plot(x, data)
      }
    }
  }

  skip("implement comparison first")


  # check the RNG will return the correct number of samples
  kumaraswamy_samples <- bayesim::rkumaraswamy(n, 0.8, 3)
  expect_equal(n, length(kumaraswamy_samples))

  # shape variable -> bound gets instable RNG, arbitrary bound instead with p_r
  test_rng(rng_fun=bayesim::rkumaraswamy, metric_mu=mean, n=n, mus=mus, shapes=ps_r,
           mu_eps=accepted_means_eps, p_acceptable_failures=p_acceptable_failures)
  # check the RNG is not too far of the input value

  # now check density function for some errors
  expect_error(bayesim::dkumaraswamy(1, 0.8)) # to few arguments
  expect_error(bayesim::dkumaraswamy(1, 0.8, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dkumaraswamy(-1, mu = 0.8, p = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dkumaraswamy("r", mu = 0.8, p = 2)) # non-numeric arguments are disallowed
  expect_error(bayesim::dkumaraswamy(1, mu = 0, p = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::dkumaraswamy(1, mu = 0.8, p = 2)) # mu is not allowed to be 1 or bigger
  expect_error(bayesim::dkumaraswamy(1, mu = 0.8, p = 0)) # p is not allowed to be 1 or smaller

  # do same for quantile function
  expect_error(bayesim::qkumaraswamy(1, 0.8)) # to few arguments
  expect_error(bayesim::qkumaraswamy(1, 0.8, 3, 4, 5)) # to many arguments
  expect_error(bayesim::qkumaraswamy("r", mu = 0.8, p = 2)) # non-numeric arguments are disallowed
  expect_error(bayesim::qkumaraswamy(1, mu = 0, p = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::qkumaraswamy(1, mu = 1, p = 2)) # mu is not allowed to be 1 or bigger
  expect_error(bayesim::qkumaraswamy(1, mu = 0.8, p = 0)) # p is not allowed to be 0 or smaller
  expect_error(bayesim::qkumaraswamy(c(-1, 2), mu = 2, p = 2)) # q is not allowed to be outside [0, 1]

  # do same for RNG function
  expect_error(bayesim::rkumaraswamy(100, 0.8)) # to few arguments
  expect_error(bayesim::rkumaraswamy(100.82, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rkumaraswamy(-1, mu = 0.8, p = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rkumaraswamy("r", mu = 0.8, p = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(bayesim::rkumaraswamy(100, mu = 0, p = 2)) # mu is not allowed to be 0 or smaller
  expect_error(bayesim::rkumaraswamy(100, mu = 1, p = 2)) # mu is not allowed to be 1 or bigger
  expect_error(bayesim::rkumaraswamy(100, mu = 0.8, p = 0)) # p is not allowed to be 0 or smaller

  expect_brms_family(n=1000, ba=0.5, int=1, shape=2, link=exp,family=bayesim::kumaraswamy,
                     rng=bayesim::rkumaraswamy, shape_name="p", thresh = 0.05)
})

