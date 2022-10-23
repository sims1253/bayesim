
n <- 100000 # number of testvalues
eps <- 1e-6
x <- seq(from = eps, to = 1 - eps, length.out = n)
x_q <- seq(from = 0.1, to = 0.9, length.out = n)

n_small <- 20
mus <- seq(from = eps, to = 1 - eps, length.out = n_small)
mus_q <- seq(from = 0.1, to = 0.9, length.out = n)
ps <- seq(from = eps, to = 20, length.out = n_small)

# mus_r <- seq(from = 1 + eps, to = 10, length.out = n_small)
ps_r <- seq(from = 0.6 + eps, to = 5, length.out = n_small)
accepted_medians_eps <- 0.02
p_acceptable_failures <- 0.05

get_q <- function(mu, p) {
  value <- -log(2) / log1p(-mu^p)
  return(value)
}

test_that("custom-kumaraswamy", {


  # calculate kumaraswamy
  dkumaraswamy_results <- dkumaraswamy(x, mu = 0.8, p = 2)
  qkumaraswamy_results <- qkumaraswamy(x, mu = 0.8, p = 2)
  # check length
  expect_equal(n, length(dkumaraswamy_results))
  expect_equal(n, length(qkumaraswamy_results))

  warning("In expect_eps quantile function, the parameters had to be quite restricted.")
  # check many shape parameters on pdf and qdf
  for (mu_idx in 1:n_small) {
    for (aux_idx in 1:n_small) {
      m <- mus[mu_idx]
      muq <- mus_q[mu_idx]
      p <- ps[aux_idx]
      pr <- ps_r[aux_idx]

      expect_eps(dkumaraswamy(x, mu = m, p = p), VGAM::dkumar(x, p, get_q(m, p)), eps, relative = TRUE)
      expect_eps(qkumaraswamy(x, mu = muq, p = pr), VGAM::qkumar(x, pr, get_q(muq, pr)), eps, relative = TRUE)
    }
  }

  # check the RNG will return the correct number of samples
  kumaraswamy_samples <- rkumaraswamy(n, 0.8, 3)
  expect_equal(n, length(kumaraswamy_samples))

  # shape variable -> bound gets instable RNG, arbitrary bound instead with p_r
  test_rng(
    rng_fun = rkumaraswamy, metric_mu = median, n = n, mu_list = mus, aux_par = ps_r,
    mu_eps = accepted_medians_eps, p_acceptable_failures = p_acceptable_failures
  )
  # check the RNG is not too far of the input value

  # now check density function for some errors
  expect_error(dkumaraswamy(1, 0.8)) # to few arguments
  expect_error(dkumaraswamy(1, 0.8, 3, 4, 5)) # to many arguments
  expect_error(dkumaraswamy(-1, mu = 0.8, p = 2)) # x is not allowed to be smaller 0
  expect_error(dkumaraswamy("r", mu = 0.8, p = 2)) # non-numeric arguments are disallowed
  expect_error(dkumaraswamy(1, mu = 0, p = 2)) # mu is not allowed to be 0 or smaller
  expect_error(dkumaraswamy(1, mu = 0.8, p = 2)) # mu is not allowed to be 1 or bigger
  expect_error(dkumaraswamy(1, mu = 0.8, p = 0)) # p is not allowed to be 1 or smaller

  # do same for quantile function
  expect_error(qkumaraswamy(1, 0.8)) # to few arguments
  expect_error(qkumaraswamy(1, 0.8, 3, 4, 5)) # to many arguments
  expect_error(qkumaraswamy("r", mu = 0.8, p = 2)) # non-numeric arguments are disallowed
  expect_error(qkumaraswamy(1, mu = 0, p = 2)) # mu is not allowed to be 0 or smaller
  expect_error(qkumaraswamy(1, mu = 1, p = 2)) # mu is not allowed to be 1 or bigger
  expect_error(qkumaraswamy(1, mu = 0.8, p = 0)) # p is not allowed to be 0 or smaller
  expect_error(qkumaraswamy(c(-1, 2), mu = 2, p = 2)) # q is not allowed to be outside [0, 1]

  # do same for RNG function
  expect_error(rkumaraswamy(100, 0.8)) # to few arguments
  expect_error(rkumaraswamy(100.82, 3, 4, 5)) # to many arguments
  expect_error(rkumaraswamy(-1, mu = 0.8, p = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(rkumaraswamy("r", mu = 0.8, p = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(rkumaraswamy(100, mu = 0, p = 2)) # mu is not allowed to be 0 or smaller
  expect_error(rkumaraswamy(100, mu = 1, p = 2)) # mu is not allowed to be 1 or bigger
  expect_error(rkumaraswamy(100, mu = 0.8, p = 0)) # p is not allowed to be 0 or smaller

  expect_error(rkuramaswamy(100, mu = 0.8, p = 1)) # kumaraswamy has to be spelled correctly!!!
  # small inside joke, given, there is a 50% chance, I misspelled it again. :P

  expect_brms_family(link = inv_logit, family = kumaraswamy, rng = rkumaraswamy, shape_name = "p")
})
