#' Title
#'
#' @param data_N
#' @param z1_x_coef
#' @param z3_x_coef
#' @param x_y_coef
#' @param z1_y_coef
#' @param z2_y_coef
#' @param y_z4_coef
#' @param x_z4_coef
#' @param y_intercept
#' @param sigma_z1
#' @param sigma_z2
#' @param sigma_z3
#' @param sigma_z4
#' @param sigma_x
#' @param sigma_y
#' @param data_link
#' @param RNG
#' @param lb
#' @param ub
#' @param oversampling
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
basedag_data <- function(data_N,
                         z1_x_coef = 1,
                         z3_x_coef = 1,
                         x_y_coef = 1,
                         z1_y_coef = 1,
                         z2_y_coef = 1,
                         y_z4_coef = 1,
                         x_z4_coef = 1,
                         y_intercept = 0,
                         sigma_z1 = 1,
                         sigma_z2 = 1,
                         sigma_z3 = 1,
                         sigma_z4 = 1,
                         sigma_x = 1,
                         sigma_y = 1,
                         data_link,
                         RNG,
                         lb,
                         ub,
                         oversampling = 10,
                         ...) {
  n <- data_N * oversampling
  z1 <- rnorm(n, sigma_z1)
  z2 <- rnorm(n, sigma_z2)
  z3 <- rnorm(n, sigma_z3)

  x <- rnorm(n, mean = (z1_x_coef * z1 + z3_x_coef * z3), sd = sigma_x)
  y <- RNG(
    n,
    data_link(
      y_intercept +
        x_y_coef * x +
        z1_y_coef * z1 +
        z2_y_coef * z2
    ),
    sigma_y
  )
  z4 <- rnorm(n, mean = (y_z4_coef * y + x_z4_coef * x), sd = sigma_z4)

  good_indices <- seq_along(y)
  if (lb) {
    good_indices <- intersect(good_indices, which(y > lb))
  }
  if (ub) {
    good_indices <- intersect(good_indices, which(y < ub))
  }
  if (length(good_indices) < data_N) {
    stop("Even with resampling, there were not enough usable samples during
         data generation. Try increasing the oversampling factor.")
  }
  good_indices <- good_indices[1:data_N]
  return(
    data.frame(
      x = x[good_indices],
      z1 = z1[good_indices],
      z2 = z2[good_indices],
      z3 = z3[good_indices],
      z4 = z4[good_indices],
      y = y[good_indices]
    )
  )
}
