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

test_that("Boolean of length n isLogic_len", {
  # correct lengths and all boolean
  expect_true(isLogic_len(TRUE))
  expect_true(isLogic_len(TRUE, 1))
  expect_true(isLogic_len(c(TRUE, FALSE), 2))
  # incorrect lengths and all boolean
  expect_false(isLogic_len(TRUE, 2))
  expect_false(isLogic_len(c(TRUE, FALSE)))
  expect_false(isLogic_len(c(TRUE, FALSE), 1))
  # correct lengths, but not boolean
  expect_false(isLogic_len("r"))
  expect_warning(expect_false(isLogic_len(isLogic_len)))
  expect_false(isLogic_len(c(TRUE, 0), 2))
  # other languages may recognise integers as boolean, which R does not do (I think)
  expect_false(isLogic_len(c("r", TRUE), 2))
  # incorrect length and non boolean
  expect_false(isLogic_len(c("r", TRUE), 3))
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
