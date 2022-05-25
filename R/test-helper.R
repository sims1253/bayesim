
#' Check that |a - b| < eps. Works with scalars and vectors on any input.
#' For vector input, one may defined, how many |a - b| >= eps are acceptable with r argument.
#'
#' @param a scalar or vector a to be compared
#' @param b scalar or vector b to be compared
#' @param eps scalar or vector, setting the max differences, eps > 0
#' @param r optional scalar (r = 0), relative number of values, that may have an difference > eps.
#' Used internally, to calculate acceptable amount of deviances. Calculated absolute value will be floored.
#'
#' @details For vector/scalar combinations, allowed are (r always scalar!):
#' \describe{
#'   \item{- all scalar}
#'   \item{- a vector, b scalar, eps scalar -> compare a values to be close enough to b}
#'   \item{- a scalar, b vector, eps scalar -> compare b values to be close enough to a}
#'   \item{- a scalar, b scalar, eps vector -> compare the difference a, b too all eps
#'   -> case not intended, but would work}
#'   \item{- a vector, b scalar, eps vector -> compare each vector-difference entry against each eps}
#'   \item{- all vectors -> each entry of |a - b| is compared to the same entry in eps}
#'   -> different vector lengths != 1 dissallowed!}
#' }
#'
#' @return success or failure with message
#'
#' @examples library(testthat)
#' library(bayesim)
#' expect_eps(1, 1.1, 0.2) # should pass
#' expect_eps(-1, -1.1, 0.2) # should pass
#' expect_eps(-1.00001, -1.00001, 1e-4) # should pass
#' expect_eps(c(0, 1, 2), c(1, 1, 2), 1e-4, 0.4) # should pass (only 1/3 were wrong, but 40% was allowed)
#' expect_eps(c(0, 1, 2), c(1, 1, 2), 1e-4, 0.1) # should fail (1/3 were wrong, but only 1/10 was allowed)
#' expect_eps(1, 1.3, 0.2) # should fail (|a - b| = 0.3 > 0.2)
#' expect_eps(c(1, 1.1), c(1, 1.1, 1.2), 0.5) # should produce an error (different vector-lengths)
#' expect_eps(1, 1.1, -0.2) # should produce an error (negative eps)
#' expect_eps(c(0, 1, 2), c(1, 1, 2), 1e-4, -0.4) # should produce an error (r too small)
#' expect_eps(c(0, 1, 2), c(1, 1, 2), 1e-4, 1.4) # should produce an error (r too big)
expect_eps <- function(a, b, eps, r = 0.0) {
  # first check, that the type is always a scalar or vector
  all_numeric <- is.numeric(a) && is.numeric(b) && is.numeric(eps) && is.numeric(r)
  if (isFALSE(all_numeric)) {
    stop("At least one of the vectors was not a scalar/vector.")
  }

  # then check, that r is only a scalar
  if (isTRUE(length(r) > 1)) {
    stop("The relative number of tolerated deviances r has to be a scalar.")
  }

  # also check, that r is in range [0, 1]
  if (isTRUE(r < 0 || r > 1)) {
    stop("The relative number of tolerated deviances r has to be in [0, 1].")
  }

  # then check, that all eps are > 0
  if (isTRUE(any(eps <= 0))) {
    stop("Tried checking against negative or zero differences. For zero difference, use expect_equal.")
  }

  # then check, if vectors are of same length
  vector_length <- max(length(a), length(b), length(eps))
  tolerated_deviances <- floor(vector_length * r)

  a_wrong <- length(a) > 1 && length(a) != vector_length
  b_wrong <- length(b) > 1 && length(b) != vector_length
  eps_wrong <- length(eps) > 1 && length(eps) != vector_length
  if (a_wrong || b_wrong || eps_wrong) {
    stop("Used different length of vectors in test.")
  }

  # last check, the actual value difference
  eps_comparison_wrong <- abs(a - b) > eps
  number_deviances <- sum(eps_comparison_wrong, na.rm = TRUE)
  if (isTRUE(number_deviances <= tolerated_deviances)) {
    succeed()
  } else {
    if (isTRUE(tolerated_deviances > 0)) {
      message <- paste("In expect_eps: ", toString(number_deviances), " of ", toString(vector_length), " were bigger, then eps. Only ", toString(tolerated_deviances), " allowed!", sep = "")
    } else {
      message <- paste("In expect_eps: ", toString(number_deviances), " of ", toString(vector_length), " were bigger, then eps. None were allowed!", sep = "")
    }
    fail(message)
  }
}

#' Check, if a is bigger than b. expect_smaller(a, b) can be done,
#' by expect_bigger(b, a) with swapped arguments.
#'
#' For vector/scalar combinations, allowed are:
#' both scalar
#' a scalar, b vector -> compare a to all b values
#' a vector, b scalar -> compare all a values to the threshold b
#' both vectors -> compare each entry to the other entry
#' -> different lengths != 1 dissallowed!
#'
#' @param a number to be tested, should be real scalar or vector
#' @param b tested against, should be real scalar or vector
#'
#' @return expect(a > threshold)
#'
#' @examples expect_bigger(1, 0) # should pass
#' expect_bigger(0, 1) # should fail
#' expect_bigger(c(1, 2), 0) # should pass
#' expect_bigger(c(0, 2), 1) # should fail
#' expect_bigger(c(1, 2), c(0, 3)) # should fail
#' expect_bigger(c(1, 2, 3), c(0, 1)) # should produce an error
expect_bigger <- function(a, b) {

  # check, that both are numerical data
  if (!is.numeric(a) || !is.numeric(b)) {
    stop("All inputs have to be numeric data.")
  }

  # check, if vectors both are vectors, both are of equal length
  vector_length <- max(length(a), length(b))
  a_wrong <- length(a) > 1 && length(a) != vector_length
  b_wrong <- length(b) > 1 && length(b) != vector_length
  if (a_wrong || b_wrong) {
    stop("Used different lengths of vectors in test")
  }

  # now compare all a and b values
  all_bigger <- all(a > b)
  if (isTRUE(all_bigger)) {
    succeed()
  } else {
    fail("At least one number a was smaller than b")
  }
}
