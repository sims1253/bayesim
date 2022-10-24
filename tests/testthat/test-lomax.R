
# Unit tests for custom lomax

# Just lambda small function, to get the alpha shape-argument, for reference lomax functions
get_lambda <- function(mu, alpha) {
  lambda <- mu * (alpha - 1)
  # lambda of extraDistr is inverted, compared to lambda in my pdf
  return(1 / lambda)
}

n <- 100000 # number of testvalues
eps <- 1e-12
x <- exp(seq(from = eps, to = 200, length.out = n)) # testset, exp(200) comes close to Max-Double
unit <- seq(from = eps, to = 1 - eps, length.out = n)

n_small <- 20
mus <- seq(from = eps, to = 1000, length.out = n_small)
alpha_list <- seq(from = 1 + eps, to = 50, length.out = n_small)

# mus_r <- seq(from = 1 + eps, to = 10, length.out = n_small)
alpha_r_list <- alpha_list + 0.2
accepted_means_eps <- 0.05
p_acceptable_failures <- 0.05

test_that("custom-lomax", {
  # calculate lomax
  dlomax_results <- dlomax(x, mu = 8, alpha = 2)
  qlomax_results <- qlomax(unit, mu = 8, alpha = 2)
  # check length
  expect_equal(n, length(dlomax_results))
  expect_equal(n, length(qlomax_results))

  # check many shape parameters on pdf and qdf
  for (m in mus) {
    for (aux_idx in seq_along(alpha_list)) {
      alpha <- alpha_list[aux_idx]
      alpha_r <- alpha_r_list[aux_idx]

      expect_eps(
        dlomax(x, mu = m, alpha = alpha),
        extraDistr::dlomax(x, get_lambda(m, alpha), alpha),
        eps,
        relative = TRUE
        )

      expect_eps(
        qlomax(unit, mu = m, alpha = alpha_r),
        extraDistr::qlomax(unit, get_lambda(m, alpha_r), alpha_r),
        eps,
        relative = TRUE
        )
    }
  }


  # check the RNG will return the correct number of samples
  lomax_samples <- rlomax(n, 1, 3)
  expect_equal(n, length(lomax_samples))

  # shape variable -> bound gets instable RNG, arbitrary bound instead with alpha_r
  test_rng(
    rng_fun = rlomax,
    metric_mu = mean,
    n = n,
    mu_list = mus,
    aux_list = alpha_r_list,
    mu_eps = accepted_means_eps,
    p_acceptable_failures = p_acceptable_failures,
    debug = TRUE,
    relative = TRUE
  )
  # check the RNG is not too far of the input value

  # now check density function for some errors
  expect_error(dlomax(1, 2)) # to few arguments
  expect_error(dlomax(1, 2, 3, 4, 5)) # to many arguments
  expect_error(dlomax(-1, mu = 2, alpha = 2)) # x is not allowed to be smaller 0
  expect_error(dlomax("r", mu = 2, alpha = 2)) # non-numeric arguments are disallowed
  expect_error(dlomax(1, mu = 0, alpha = 2)) # mu is not allowed to be 0 or smaller
  expect_error(dlomax(1, mu = 1, alpha = 0)) # alpha is not allowed to be 1 or smaller

  # do same for quantile function
  expect_error(qlomax(1, 2)) # to few arguments
  expect_error(qlomax(1, 2, 3, 4, 5)) # to many arguments
  expect_error(qlomax("r", mu = 2, alpha = 2)) # non-numeric arguments are disallowed
  expect_error(qlomax(1, mu = 0, alpha = 2)) # mu is not allowed to be 0 or smaller
  expect_error(qlomax(1, mu = 1, alpha = 0)) # alpha is not allowed to be 0 or smaller
  expect_error(qlomax(c(-1, 2), mu = 2, alpha = 2)) # q is not allowed to be outside [0, 1]

  # do same for RNG function
  expect_error(rlomax(100, 2)) # to few arguments
  expect_error(rlomax(10, 2, 3, 4, 5)) # to many arguments
  expect_error(rlomax(-1, mu = 2, alpha = 2)) # number of drawn samples cannot be smaller 0
  expect_warning(expect_error(rlomax("r", mu = 2, alpha = 2))) # non-numeric arguments are disallowed
  # also non-numeric arguments for n will throw warning
  expect_error(rlomax(100, mu = 0, alpha = 2)) # mu is not allowed to be 0 or smaller
  expect_error(rlomax(100, mu = 1, alpha = 0)) # alpha is not allowed to be 0 or smaller

  expect_brms_family(
    intercept = 5,
    ref_intercept = 5,
    aux_par = 2,
    link = exp,
    family = lomax,
    rng = rlomax,
    aux_name = "alpha"
    )
})
