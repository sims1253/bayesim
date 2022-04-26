library(testthat)
#' Title
#'
#' @param a scalar or vector a to be compared
#' @param b scalar or vector b to be compared
#' @param eps scalar or vector, setting the max difference
#'
#' @return expect_true(a close enought to b)
#' @export
#'
#' @examples expect_eps(1, 1.1, 0.2) should pass
#' expect_eps(1, 1.3, 0.2) should fail
#' a <- {1, 2}
#' b <- {3, 4, 5}
#' expect_eps(a, b, 10) should produce an error
expect_eps <- function(a, b, eps) {
  # first check, if vectors are of same length
  vector_length <- max(length(a), length(b), length(eps))

  a_wrong   <- length(a) > 1 && length(a) != vector_length
  b_wrong   <- length(b) > 1 && length(b) != vector_length
  eps_wrong <- length(eps) > 1 && length(eps) != vector_length
  if(a_wrong || b_wrong || eps_wrong)
    stop("Used different length of vectors in test")

  bool_comparison <- abs(a - b) < eps
  if(isTRUE(all(bool_comparison)))
    succeed("All entries of a and be were in the expected range")
  else
    fail("At least one comparison between a and b had a bigger difference, then eps")
}
