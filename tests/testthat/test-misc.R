
eps <- 1e-6
n_testset <- 100

# Test the link functions
test_that("logit-link", {
  # check, that the error is within reason
  expect_eps(0, logit(0.5), eps)
  expect_eps(-1.38629436112, logit(0.2), eps)
  expect_eps(0.847297860387, logit(0.7), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(0, -1.38629436112, 0.847297860387), logit(c(0.5, 0.2, 0.7)), eps)
  # check link against link from another package
  testset_logit <- seq(from = eps, to = 1 - eps, length.out = n_testset)
  expect_eps(brms:::logit(testset_logit), logit(testset_logit), eps)
  # check values at the boundary of the defined space
  expect_equal(logit(0), -Inf)
  expect_equal(logit(1), Inf)
  # check, that non-numeric arguments result in an error
  expect_error(logit("R"))
  # check, that wrong number of arguments produce an error
  expect_error(logit(0.1, 0.5))
  # check for NaN and warning
  expect_warning(expect_true(is.na(logit(-1))))
  expect_warning(expect_true(is.na(logit(2))))
})

test_that("inverse-logit-link", {
  # check, that the error is within reason
  expect_eps(0.5, inv_logit(0), eps)
  expect_eps(0.119202922022, inv_logit(-2), eps)
  expect_eps(0.73105857863, inv_logit(1), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(0.5, 0.119202922022, 0.73105857863), inv_logit(c(0, -2, 1)), eps)
  # check link against link from another package
  testset_invlogit <- seq(from = -10, to = 10, length.out = n_testset)
  expect_eps(brms:::inv_logit(testset_invlogit), inv_logit(testset_invlogit), eps)
  # check values, for x to Inf
  expect_equal(1, inv_logit(Inf))
  expect_equal(0, inv_logit(-Inf))
  # check, that non-numeric arguments result in an error
  expect_error(inv_logit("R"))
  # check, that wrong number of arguments produce an error
  expect_error(logit(0.1, 0.5))
})

test_that("logistic-link", {
  # same tests as inverse-logit-link, given, it another name for the same link
  # check, that the error is within reason
  expect_eps(0.5, logistic(0), eps)
  expect_eps(0.119202922022, logistic(-2), eps)
  expect_eps(0.73105857863, logistic(1), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(0.5, 0.119202922022, 0.73105857863), logistic(c(0, -2, 1)), eps)
  # equivalent to inverse logit
  testset_logistic <- seq(from = -100, to = 100, length.out = n_testset)
  expect_equal(inv_logit(testset_logistic), logistic(testset_logistic))
  # also check against other package
  expect_eps(brms:::inv_logit(testset_logistic), logistic(testset_logistic), eps)
  # check values, for x to Inf/x on boundary of defined space
  expect_equal(1, logistic(Inf))
  expect_equal(0, logistic(-Inf))
  # check, that non-numeric arguments result in an error
  expect_error(logistic("R"))
  # check, that wrong number of arguments produce an error
  expect_error(logistic(0.1, 0.5))
})

test_that("cloglog-link", {
  # check, that the error is within reason
  expect_eps(-9.21029, cloglog(0.0001), eps)
  expect_eps(-2.250367, cloglog(0.1), eps)
  expect_eps(-0.005764308, cloglog(0.63), eps)
  expect_eps(0.475885, cloglog(0.8), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(-9.21029, -2.250367, -0.005764308, 0.475885), cloglog(c(0.0001, 0.1, 0.63, 0.8)), eps)
  # check link against link from another package
  testset_cloglog <- seq(from = eps, to = 1 - eps, length.out = n_testset)
  expect_eps(brms:::cloglog(testset_cloglog), cloglog(testset_cloglog), eps)
  # check boundary values
  expect_equal(-Inf, cloglog(0))
  expect_equal(Inf, cloglog(1))
  # check, that non-numeric arguments result in an error
  expect_error(cloglog("R"))
  # check, that wrong number of arguments produce an error
  expect_error(cloglog(0.1, 0.5))
  # check values outside the defined scope are NaN and throw warning
  expect_warning(expect_true(is.na(cloglog(-1))))
  expect_warning(expect_true(is.na(cloglog(2))))
})

test_that("inverse-cloglog-link", {
  # check, that the error is within reason
  expect_eps(0.04856801, inv_cloglog(-3), eps)
  expect_eps(0.3077994, inv_cloglog(-1), eps)
  expect_eps(0.6321206, inv_cloglog(0), eps)
  expect_eps(0.934012, inv_cloglog(1), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(0.04856801, 0.3077994, 0.6321206, 0.934012), inv_cloglog(c(-3, -1, 0, 1)), eps)
  # check against another package
  testset_invcloglog <- seq(from = -10, to = 5, length.out = n_testset)
  expect_eps(brms:::inv_cloglog(testset_invcloglog), inv_cloglog(testset_invcloglog), eps)
  # check, values of x to Inf, on the boundary of the defined space
  expect_equal(0, inv_cloglog(-Inf))
  expect_equal(1, inv_cloglog(Inf))
  # check, that non-numeric arguments result in an error
  expect_error(inv_cloglog("R"))
  # check, that wrong number of arguments produce an error
  expect_error(inv_cloglog(0.1, 0.5))
})

test_that("cauchit-link", {
  # check, that the error is within reason
  expect_eps(-3.077684, cauchit(0.1), eps)
  expect_eps(-0.3249197, cauchit(0.4), eps)
  expect_eps(0, cauchit(0.5), eps)
  expect_eps(1.376382, cauchit(0.8), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(-3.077684, -0.3249197, 0, 1.376382), cauchit(c(0.1, 0.4, 0.5, 0.8)), eps)
  # check against another package
  testset_cauchit <- seq(from = eps, to = 1 - eps, length.out = n_testset)
  expect_eps(qcauchy(testset_cauchit), cauchit(testset_cauchit), 10 * eps)
  # check values on the boundary, should logically be Inf, but just approach Inf
  expect_equal(cauchit(0), -Inf)
  expect_equal(cauchit(1), Inf)
  # check values outside the defined scope are NaN and throw warning
  expect_warning(expect_true(is.na(cauchit(-1))))
  expect_warning(expect_true(is.na(cauchit(2))))
  # check, that non-numeric arguments result in an error
  expect_error(cauchit("R"))
  # check, that wrong number of arguments produce an error
  expect_error(cauchit(0.1, 0.5))
})

test_that("inv-cauchit-link", {
  # check, that the error is within reason
  expect_eps(0.03172552, inv_cauchit(-10), eps)
  expect_eps(0.5, inv_cauchit(0), eps)
  expect_eps(0.9682745, inv_cauchit(10), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(0.03172552, 0.5, 0.9682745), inv_cauchit(c(-10, 0, 10)), eps)
  # equivalent to pcauchy
  testset_invcauchit <- seq(from = -10, to = 10, length.out = n_testset)
  expect_equal(pcauchy(testset_invcauchit), inv_cauchit(testset_invcauchit))
  # check x approaching Inf on boundry of defined space
  expect_equal(0, inv_cauchit(-Inf))
  expect_equal(1, inv_cauchit(Inf))
  # check, that non-numeric arguments result in an error
  expect_error(inv_cauchit("R"))
  # check, that wrong number of arguments produce an error
  expect_error(inv_cauchit(0.1, 0.5))
})

test_that("gaussian-error-function", {
  # check, that the error is within reason
  expect_eps(-0.8427008, erf(-1), eps)
  expect_equal(0, erf(0))
  expect_eps(0.2227026, erf(0.2), eps)
  # check vector as argument returns vector with same results
  expect_eps(c(-0.8427008, 0, 0.2227026), erf(c(-1, 0, 0.2)), eps)
  # check x approaching Inf on boundry of defined space
  expect_equal(-1, erf(-Inf))
  expect_equal(1, erf(Inf))
  # check, that non-numeric arguments result in an error
  expect_error(erf("R"))
  # check, that wrong number of arguments produce an error
  expect_error(erf(0.1, 0.5))
})

test_that("Softplus link-function", {
  # check, that the error is within reason
  expect_eps(-2.252168, softplus(0.1), eps)
  expect_eps(0.5413249, softplus(1), eps)
  expect_eps(3.981515, softplus(4), eps)
  expect_equal(100, softplus(100)) # should be equal to machine precision (I think)
  # check vector as argument returns vector with same results
  expect_eps(c(-2.252168, 0.5413249, 3.981515, 100), softplus(c(0.1, 1, 4, 100)), eps)
  # check x approaching Inf on boundry of defined space
  expect_equal(-Inf, softplus(0))
  expect_equal(Inf, softplus(Inf))
  # check, that non-numeric arguments result in an error
  expect_error(softplus("R"))
  # check, that wrong number of arguments produce an error
  expect_error(softplus(0.1, 0.5))
  # check values outside the defined scope are NaN and throw warning
  expect_warning(expect_true(is.na(softplus(-1))))
})

test_that("Inverse Softplus link-function", {
  # check, that the error is within reason
  expect_eps(4.53989e-05, inv_softplus(-10), eps)
  expect_eps(0.3132617, inv_softplus(-1), eps)
  expect_eps(1.313262, inv_softplus(1), eps)
  expect_equal(100, inv_softplus(100)) # should be equal to machine precision (I think)
  # check vector as argument returns vector with same results
  expect_eps(c(4.53989e-05, 0.3132617, 1.313262, 100), inv_softplus(c(-10, -1, 1, 100)), eps)
  # check x approaching Inf on boundry of defined space
  expect_equal(0, inv_softplus(-Inf))
  expect_equal(Inf, inv_softplus(Inf))
  # check, that non-numeric arguments result in an error
  expect_error(inv_softplus("R"))
  # check, that wrong number of arguments produce an error
  expect_error(inv_softplus(0.1, 0.5))
})

# Test a few fúrther miscellanios helping functions
test_that("Real of length n isNum_len", {
  # correct lengths and all numerics
  expect_true(isNum_len(0.2))
  expect_true(isNum_len(0.2, 1))
  expect_true(isNum_len(c(0.2, 0.3), 2))
  # incorrect lengths and all numerics
  expect_false(isNum_len(0.2, 2))
  expect_false(isNum_len(c(0.2, 0.3)))
  expect_false(isNum_len(c(0.2, 0.3), 1))
  # correct lengths, but not numerics
  expect_false(isNum_len("r"))
  expect_warning(expect_false(isNum_len(isNum_len)))
  expect_false(isNum_len(c(0.2, NA), 2)) # numeric vectors containing NAs are
  # usually also called numeric, but should not be set as such
  expect_false(isNum_len(c("r", 0.2), 2))
  expect_false(isNum_len(c("r", 0.2), 2))
  # incorrect length and non numerics
  expect_false(isNum_len(c("r", 0.2), 3))
})

test_that("Integer of length n isInt_len", {
  # correct lengths and all integers
  expect_true(isInt_len(1))
  expect_true(isInt_len(1, 1))
  expect_true(isInt_len(c(1, 2), 2))
  # incorrect lengths and all integers
  expect_false(isInt_len(1, 2))
  expect_false(isInt_len(c(1, 2)))
  expect_false(isInt_len(c(1, 2), 1))
  # correct lengths, but not integers
  expect_false(isInt_len("r"))
  expect_false(isInt_len(c("r", 1), 2))
  expect_false(isInt_len(c(1.1, 1), 2))
  expect_warning(expect_false(isInt_len(isInt_len)))
  expect_false(isInt_len(c(1, NA), 2)) # numeric vectors containing NAs are
  # usually also called numeric, but should not be set as such
  expect_false(isInt_len(c("r", 1), 2))
  # incorrect lengths and non integers
  expect_false(isInt_len(c("r", 2), 3))
})

test_that("Single string isSingleString", {
  # correct strings
  expect_true(isSingleString("r"))
  expect_true(isSingleString("abc"))
  expect_true(isSingleString(c("abc")))
  # correct strings, but not single string
  expect_false(isSingleString(c("abc", "r")))
  # not a string
  expect_false(isSingleString(1))
  expect_warning(expect_false(isSingleString(isSingleString)))
  expect_false(isSingleString(c(NA, "abc")))
  # multiple inputs containing non strings
  expect_false(isSingleString(c("abc", 1, isSingleString)))
})

test_that("data limiting function limit_data", {
  input <- c(1, 2, 3, 4)
  result_both <- limit_data(input, c(2, 3))
  # length should not be changed
  expect_true(length(result_both) == 4)
  expect_true(all(2 <= result_both | result_both <= 3))
  # now only change the upper boundary
  result_upper <- limit_data(input, c(NA, 3))
  expect_true(all(1 <= result_upper | result_upper <= 3))
  # also the first input should not be changed now
  expect_true(input[1] == result_upper[1])
  # ... while the last one was set to 3
  expect_true(3 == result_upper[4])
  # should not work for limits beeing wrong way around
  expect_error(limit_data(input, c(3, 2)))
  # wrong amount of limit arguments
  expect_error(limit_data(input, c(1, 2, 3)))
  expect_error(limit_data(input, c(1)))
  # wrong input data
  expect_error(limit_data(c(1, "r"), c(1, 2)))
  expect_error(limit_data(c(1, NA), c(1, 2)))
})
