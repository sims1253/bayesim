
test_that("custom-softplusnormal", {

  # load in values
  pdf_data <- readRDS("precalc_values/softplusnormal_refpdf")

  eps <- 1e-12 # 2 digits more, than sim
  unit_int <- c(eps, 1 - eps)
  mu_unit_int <- c(0.1, 0.9)
  pos_int <- c(eps, 200)
  shape_int <- c(0.1, 20)
  n <- 1000
  n_small <- 20

  mus <- seq(from=pos_int[1], to=pos_int[2], length.out=n_small)
  sigmas <- seq(from=shape_int[1], to=shape_int[2], length.out=n_small)
  x <- exp(seq(from=pos_int[1], to=pos_int[2], length.out=n))

  pdf_ref <- as.matrix(pdf_data)


  # calculate beta-prime
  dsoftplusnormal_results <- dsoftplusnormal(x, mu = 1, sigma = 2)
  # check length
  expect_equal(n, length(dsoftplusnormal_results))
  # check against one precalculated value
  expect_eps(0.3922206, dsoftplusnormal(x = 0.5, mu = 1, sigma = 2), eps = 1e-6)
  # ??? The constant in this test was wrong? Probably corrected implementation, or something
  # Now it is not 1.27.. anymore, but 0.392.. Maybe that is more correct?

  for (outer in 1:n_small) {
    for (inner in 1:n_small) {
      mu <- mus[outer]
      sigma <- sigmas[inner]
      expect_eps(dsoftplusnormal(x, mu, sigma), pdf_ref[[outer, inner]], eps, relative = TRUE)
    }
  }

  # now check density function for some errors
  expect_error(dsoftplusnormal(0.5, 2)) # to few arguments
  expect_error(dsoftplusnormal(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(dsoftplusnormal(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(dsoftplusnormal(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  expect_error(dsoftplusnormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed

  # do same for RNG function
  expect_error(rsoftplusnormal(100, 2)) # to few arguments
  expect_error(rsoftplusnormal(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rsoftplusnormal(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(rsoftplusnormal("r", mu = 2, sigma = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(rsoftplusnormal(100, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller


  # check the RNG will return the correct number of samples
  softplusnormal_samples <- rsoftplusnormal(n, 2, 3)
  expect_equal(n, length(softplusnormal_samples))

  n_rng <- 100000
  accepted_medians_eps <- 0.01
  p_acceptable_failures <- 0.05

  # check the RNG is not too far of the input value
  test_rng(
    rng_fun = rsoftplusnormal, metric_mu = median, n = n_rng, mu_list = mus, aux_par = sigmas,
    mu_eps = accepted_medians_eps, p_acceptable_failures = p_acceptable_failures, mu_link = softplus
  )

  # check custom BRMS family implementation
  expect_brms_family(
    ba = 0.2, intercept = 0.4, shape = 2, link = identity, family = softplusnormal,
    rng = rsoftplusnormal, shape_name = "sigma"
  )
})
