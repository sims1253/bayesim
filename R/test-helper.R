
#' Check that |a - b| < eps. Works with scalars and vectors on any input.
#'
#' @param a scalar or vector a to be compared
#' @param b scalar or vector b to be compared
#' @param eps scalar or vector, setting the max differences, eps > 0
#'
#' @return expect_true(a close enought to b)
#' @export
#'
#' @examples expect_eps(1, 1.1, 0.2) should pass
#' expect_eps(-1, -1.1, 0.2) should pass
#' expect_eps(-1.00001, -1.00001, 1e-4) should pass
#' expect_eps(1, 1.3, 0.2) should fail
#' expect_eps(c(1, 1.1), c(1, 1.1, 1.2), 0.5) should produce an error
#' expect_eps(1, 1.1, -0.2) should produce an error
expect_eps <- function(a, b, eps) {
  # first check, that the type is always a scalar or vector
  all_numeric <- is.numeric(a) && is.numeric(b) && is.numeric(eps)
  if(isFALSE(all_numeric))
    stop("At least one of the vectors was not a scalar/vector.")

  # then check, that
  # first check, that all eps are > 0
  if(isTRUE(any(eps <= 0)))
    stop("Tried checking against negative/zero differences")

  # then check, if vectors are of same length
  vector_length <- max(length(a), length(b), length(eps))

  a_wrong   <- length(a) > 1 && length(a) != vector_length
  b_wrong   <- length(b) > 1 && length(b) != vector_length
  eps_wrong <- length(eps) > 1 && length(eps) != vector_length
  if(a_wrong || b_wrong || eps_wrong)
    stop("Used different length of vectors in test")

  # last check, the actual value difference
  bool_comparison <- abs(a - b) < eps
  if(isTRUE(all(bool_comparison)))
    succeed()
  else
    fail("At least one comparison between a and b had a bigger difference, then eps")
}

#' Check, if a is bigger than b. expect_smaller(a, b) can be done,
#' by expect_bigger(b, a) with swapped arguments.
#'
#' @param a number to be tested, should be real scalar or vector
#' @param b tested against, should be real scalar or vector
#'
#' @return expect(a > threshold)
#' @export
#'
#' @examples expect_bigger(1, 0) should pass
#' expect_bigger(0, 1) should fail
#' expect_bigger(c(1, 2), 0) should pass
#' expect_bigger(c(0, 2), 1) should fail
#' expect_bigger(c(1, 2), c(0, 3)) should fail
#' expect_bigger(c(1, 2, 3), c(0, 1)) should produce an error
expect_bigger <- function(a, b) {
  if(!is.numeric(a) || !is.numeric(b))
    stop("All inputs have to be numeric data.")

  vector_length <- max(length(a), length(b))
  a_wrong <- length(a) > 1 && length(a) != vector_length
  b_wrong <- length(b) > 1 && length(b) != vector_length
  if(a_wrong || b_wrong)
    stop("Used different length of vectors in test")

  bool_comparison <- all(a > b)
  if(isTRUE(bool_comparison))
    succeed()
  else
    fail("At least one number was smaller, than b")

}
