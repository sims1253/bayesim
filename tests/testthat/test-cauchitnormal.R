
test_that("custom-cauchitnormal", {

  # load in values
  data <- readRDS("precalc_values/cauchitnormal_refdata")
  pdf_data <- readRDS("precalc_values/cauchitnormal_refpdf")

  n <- data$n
  n_small <- data$n_small
  eps <- data$eps
  mus <- data$mus
  sigmas <- data$shapes

  x <- data$x
  pdf_ref <- as.matrix(pdf_data)

  # calculate beta-prime
  dcauchitnormal_results <- bayesim::dcauchitnormal(x, mu = 1, sigma = 2)
  # check length
  expect_equal(n, length(dcauchitnormal_results))
  # check against one precalculated value
  expect_eps(0.5530229, bayesim::dcauchitnormal(x=0.5, mu=1, sigma=2), eps)

  # check the RNG will return the correct number of samples
  cauchitnormal_samples <- bayesim::rcauchitnormal(n, 2, 3)
  expect_equal(n, length(cauchitnormal_samples))

  for(outer in 1:n_small) {
    for(inner in 1:n_small) {
      mu <- mus[outer]
      sigma <- sigmas[inner]
      expect_eps(bayesim::dcauchitnormal(x, mu, sigma), pdf_ref[[outer, inner]], eps)
    }
  }

  accepted_medians_eps <- 0.15
  p_acceptable_failures <- 0.05
  n_rng <- 100000

  # check the RNG is not too far of the input value p_acceptable_failures
  test_rng(rng_fun=bayesim::rcauchitnormal, metric_mu=median, n=n_rng, mus=mus, shapes=sigmas,
           mu_eps=accepted_medians_eps, p_acceptable_failures=p_acceptable_failures, mu_link=cauchit)



  # now check density function for some errors
  expect_error(bayesim::dcauchitnormal(0.5, 2)) # to few arguments
  expect_error(bayesim::dcauchitnormal(0.5, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::dcauchitnormal(-1, mu = 2, sigma = 2)) # x is not allowed to be smaller 0
  expect_error(bayesim::dcauchitnormal(2, mu = 2, sigma = 2)) # x is not allowed to be bigger 1
  expect_error(bayesim::dcauchitnormal(0.5, mu = 1, sigma = -1)) # sigma is not allowed to be 0 or smaller
  expect_error(bayesim::dcauchitnormal("r", mu = 2, sigma = 2)) # non-numeric arguments are disallowed


  # do same for RNG function
  expect_error(bayesim::rcauchitnormal(100, 2)) # to few arguments
  expect_error(bayesim::rcauchitnormal(10, 2, 3, 4, 5)) # to many arguments
  expect_error(bayesim::rcauchitnormal(-1, mu = 2, sigma = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(bayesim::rcauchitnormal("r", mu = 2, sigma = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning

  expect_brms_family(ba=0.2, intercept=0.4, shape=2, link=identity, family=bayesim::cauchitnormal,
                     rng=bayesim::rcauchitnormal, shape_name="sigma")
})
