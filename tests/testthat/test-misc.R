
eps <- 1e-6

test_that("logit-link", {
  # check, that the error is within reason
  expect_eps(0, logit(0.5), eps)
  expect_eps(-1.38629436112, logit(0.2), eps)
  expect_eps(0.847297860387, logit(0.7), eps)
  # check values outside the defined space, to throw warning
  expect_warning(logit(-1))
  expect_warning(logit(2))
  # check values at the boundary of the defined space
  expect_equal(logit(0), -Inf)
  expect_equal(logit(1), Inf)
  # check, that non-numeric arguments result in an error
  expect_error(logit("R"))
  expect_error(logit("R", 0.5))
})

test_that("inverse-logit-link", {
  # check, that the error is within reason
  expect_eps(0.5, inv_logit(0), eps)
  expect_eps(0.119202922022, inv_logit(-2), eps)
  expect_eps(0.73105857863, inv_logit(1), eps)
  # check values, for x to Inf
  expect_eps(1, inv_logit(Inf), eps)
  expect_eps(0, inv_logit(-Inf), eps)
  # check, that non-numeric arguments result in an error
  expect_error(inv_logit("R"))
  expect_error(inv_logit("R", 0.5))
})

test_that("logistic-link", {
  # same tests as inverse-logit-link, given, it another name for the same link
  # check, that the error is within reason
  expect_eps(0.5, logistic(0), eps)
  expect_eps(0.119202922022, logistic(-2), eps)
  expect_eps(0.73105857863, logistic(1), eps)
  # check values, for x to Inf/x on boundary of defined space
  expect_eps(1, logistic(Inf), eps)
  expect_eps(0, logistic(-Inf), eps)
  # check, that non-numeric arguments result in an error
  expect_error(logistic("R"))
  expect_error(logistic("R", 0.5))
})

test_that("cloglog-link", {
  # check, that the error is within reason
  expect_eps(-9.21029, cloglog(0.0001), eps)
  expect_eps(-2.250367, cloglog(0.1), eps)
  expect_eps(-0.005764308, cloglog(0.63), eps)
  expect_eps(0.475885, cloglog(0.8), eps)
  # check boundary values
  expect_equal(-Inf, cloglog(0))
  expect_equal(Inf, cloglog(1))
  # check values outside of the defined space to throw warning
  expect_warning(logit(-1))
  expect_warning(logit(2))
  # check, that non-numeric arguments result in an error
  expect_error(cloglog("R"))
  expect_error(cloglog("R", 0.5))
})

test_that("inverse-cloglog-link", {
  # check, that the error is within reason
  expect_eps(0.04856801, inv_cloglog(-3), eps)
  expect_eps(0.3077994, inv_cloglog(-1), eps)
  expect_eps(0.6321206, inv_cloglog(0), eps)
  expect_eps(0.934012, inv_cloglog(1), eps)
  # check, values of x to Inf, on the boundary of the defined space
  expect_equal(0, inv_cloglog(-Inf))
  expect_equal(1, inv_cloglog(Inf))
  # check, that non-numeric arguments result in an error
  expect_error(inv_cloglog("R"))
  expect_error(inv_cloglog("R", 0.5))
})

test_that("cauchit-link", {
  # check, that the error is within reason
  expect_eps(-3.077684, cauchit(0.1), eps)
  expect_eps(-0.3249197, cauchit(0.4), eps)
  expect_eps(0, cauchit(0.5), eps)
  expect_eps(1.376382, cauchit(0.8), eps)
  # cauchit is periodic, so check that
  expect_eps(cauchit(1.1), cauchit(0.1), eps)
  expect_eps(cauchit(1.4), cauchit(0.4), eps)
  expect_eps(cauchit(1.5), cauchit(0.5), eps)
  expect_eps(cauchit(1.8), cauchit(0.8), eps)
  # check values on the boundary, should logically be Inf, but just approach Inf
  expect_bigger(abs(cauchit(0)), 1e+15)
  expect_bigger(abs(cauchit(1)), 1e+15)
  # check, that non-numeric arguments result in an error
  expect_error(inv_cloglog("R"))
  expect_error(inv_cloglog("R", 0.5))
})

test_that("cauchit-link", {
  # check, that the error is within reason
  expect_eps(0.03172552, inv_cauchit(-10), eps)
  expect_eps(0.5, inv_cauchit(0), eps)
  expect_eps(0.9682745, inv_cauchit(10), eps)
  # check x approaching Inf on boundry of defined space
  expect_eps(0, inv_cauchit(-Inf), eps)
  expect_eps(1, inv_cauchit(Inf), eps)
  # check, that non-numeric arguments result in an error
  expect_error(inv_cauchit("R"))
  expect_error(inv_cauchit("R", 0.5))
})
