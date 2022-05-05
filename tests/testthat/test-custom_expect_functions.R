
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
  expect_success(expect_eps(c(2, 1, 1), 1.1, eps=0.2, r=0.4))
  expect_failure(expect_eps(c(1, 2), 1.1, 0.2))
  expect_failure(expect_eps(c(1, 1.1), 2, 0.2))
  expect_failure(expect_eps(c(2, 1, 1), 1.1, eps=0.2, r=0.2))
  # a vector, b vector, eps scalar
  expect_success(expect_eps(c(1, 1), c(1.1, 1.1), 0.2))
  expect_success(expect_eps(c(2, 1, 1), c(1.1, 1.1, 1.1), eps=0.2, r=0.4))
  expect_failure(expect_eps(c(1, 2), c(1.1, 1.1), 0.2))
  expect_failure(expect_eps(c(2, 1, 1), c(1.1, 1.1, 1.1), eps=0.2, r=0.2))
  expect_error(expect_eps(c(1, 1, 1), c(1.1, 1.1), 0.2))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1), 0.2))
  # a vector, b scalar, eps vector
  expect_success(expect_eps(c(1, 1), 1.1, c(0.2, 0.3)))
  expect_success(expect_eps(c(1, 1, 2), 1.1, eps=c(0.2, 0.3, 0.2), r=0.4))
  expect_success(expect_eps(c(1, 1), 1.1, c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 2), 1.1, c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 1), 2, c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 1, 2), 1.1, eps=c(0.2, 0.3, 0.2), r=0.2))
  expect_error(expect_eps(c(1, 1, 1), 1.1, c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 2, 3), 1.1, c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 1), 1.1, c(0.2, -0.3)))
  # all vectors
  expect_success(expect_eps(c(1, 1), c(1.1, 1.1), c(0.2, 0.3)))
  expect_success(expect_eps(c(1, 3), c(1.1, 3.1), c(0.2, 0.3)))
  expect_success(expect_eps(c(1, 5), c(1.1, 7), c(0.2, 3)))
  expect_success(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps=c(0.2, 0.3, 0.2), r=0.4))
  expect_failure(expect_eps(c(1, 2), c(1.1, 1.1), c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps=c(0.2, 0.3, 0.2), r=0.2))
  expect_error(expect_eps(c(1, 1, 1), c(1.1, 1.1), c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1), c(0.2, 0.3)))
  # additional r-error-tests
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps=c(0.2, 0.3, 0.2), r=-0.1))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps=c(0.2, 0.3, 0.2), r=1.1))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1, 1.1), eps=c(0.2, 0.3, 0.2), r=c(0.1, 0.2, 0.3)))
  expect_error(expect_eps(1, 1.1, 0.2, "r"))
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
