test_that("test custom expect_eps function", {
  # wrong amount of arguments
  expect_error(expect_eps(1, 1))
  expect_error(expect_eps(1, 1, 1, 1))
  # all scalars
  expect_success(expect_eps(1, 1.1, 0.2))
  expect_success(expect_eps(1, 3, 3))
  expect_success(expect_eps(-1, -1.1, 0.2))
  expect_failure(expect_eps(1, 2, 0.2))
  expect_failure(expect_eps(2, 1, 0.2))
  expect_error(expect_eps(0.1, 0.11, -0.1))
  # a vector, b scalar, eps scalar
  expect_success(expect_eps(c(1, 1), 1.1, 0.2))
  expect_failure(expect_eps(c(1, 2), 1.1, 0.2))
  expect_failure(expect_eps(c(1, 1.1), 2, 0.2))
  # a vector, b vector, eps scalar
  expect_success(expect_eps(c(1, 1), c(1.1, 1.1), 0.2))
  expect_failure(expect_eps(c(1, 2), c(1.1, 1.1), 0.2))
  expect_error(expect_eps(c(1, 1, 1), c(1.1, 1.1), 0.2))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1), 0.2))
  # a vector, b scalar, eps vector
  expect_success(expect_eps(c(1, 1), 1.1, c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 2), 1.1, c(0.2, 0.3)))
  expect_failure(expect_eps(c(1, 1), 2, c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 1, 1), 1.1, c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 2, 3), 1.1, c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 1), 1.1, c(0.2, -0.3)))
  # all vectors
  expect_success(expect_eps(c(1, 1), c(1.1, 1.1), c(0.2, 0.3)))
  expect_success(expect_eps(c(1, 3), c(1.1, 3.1), c(0.2, 0.3)))
  expect_success(expect_eps(c(1, 5), c(1.1, 7), c(0.2, 3)))
  expect_failure(expect_eps(c(1, 2), c(1.1, 1.1), c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 1, 1), c(1.1, 1.1), c(0.2, 0.3)))
  expect_error(expect_eps(c(1, 1, 2), c(1.1, 1.1), c(0.2, 0.3)))
})
