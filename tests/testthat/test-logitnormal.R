
test_that("custom-logitnormal", {

  # load in values
  pdf_data <- readRDS("precalc_values/logitnormal_refpdf")

  eps <- 1e-12 # 2 digits more, than sim
  unit_int <- c(eps, 1 - eps)
  mu_unit_int <- c(0.1, 0.9)
  pos_int <- c(eps, 200)
  shape_int <- c(0.1, 20)
  n <- 1000
  n_small <- 20

  mus <- logit(seq(from=mu_unit_int[1], to=mu_unit_int[2], length.out=n_small))
  sigmas <- seq(from=shape_int[1], to=shape_int[2], length.out=n_small)
  x <- seq(from=unit_int[1], to=unit_int[2], length.out=n)

  pdf_ref <- as.matrix(pdf_data)

  # calculate beta-prime
  dlogitnormal_results <- bayesim::dlogitnormal(x, mu = 1, sigma = 2)
  # check length
  expect_equal(n, length(dlogitnormal_results))
  # check against one precalculated value
  expect_eps(0.7041307, bayesim::dlogitnormal(x = 0.5, mu = 1, sigma = 2), 1e-6)

  # check the RNG will return the correct number of samples
  logitnormal_samples <- bayesim::rlogitnormal(n, 2, 3)
  expect_equal(n, length(logitnormal_samples))

  for (outer in 1:n_small) {
    for (inner in 1:n_small) {
      mu <- mus[outer]
      sigma <- sigmas[inner]
      expect_eps(bayesim::dlogitnormal(x, mu, sigma), pdf_ref[[outer, inner]], eps, relative = TRUE)
    }
  }

  n_rng <- 100000
  accepted_medians_eps <- 0.2
  p_acceptable_failures <- 0.1
  # check the RNG is not too far of the input value
  test_rng(
    rng_fun = bayesim::rlogitnormal, metric_mu = median, n = n_rng, mu_list = mus, aux_par = sigmas,
    mu_eps = accepted_medians_eps, p_acceptable_failures = p_acceptable_failures, mu_link = logit
  )



  # now check density function for some errors
  expect_error(bayesim::dlogitnormal(0.5, 2)) # to few arguments
  expect_error(bayesim::dlogitnormal(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dlogitnormal(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dlogitnormal(2, mu = 2, sigma = 2)) # x is not allowed to be bigger 1
  expect_error(bayesim::dlogitnormal(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  expect_error(bayesim::dlogitnormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed


  # do same for RNG function
  expect_error(bayesim::rlogitnormal(100, 2)) # to few arguments
  expect_error(bayesim::rlogitnormal(10, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rlogitnormal(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rlogitnormal("r", mu = 2, sigma = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(bayesim::rlogitnormal(100, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller

  expect_brms_family(
    ba = 0.2, intercept = 0.4, shape = 2, link = identity, family = bayesim::logitnormal,
    rng = bayesim::rlogitnormal, shape_name = "sigma"
  )
})
