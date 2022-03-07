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
#' @param lb
#' @param ub
#' @param ...
#' @param data_family
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
basedag_data <- function(data_N,
                         z1_x_coef,
                         z3_x_coef,
                         x_y_coef,
                         z1_y_coef,
                         z2_y_coef,
                         y_z4_coef,
                         x_z4_coef,
                         y_intercept,
                         sigma_z1,
                         sigma_z2,
                         sigma_z3,
                         sigma_z4,
                         sigma_x,
                         sigma_y,
                         data_link,
                         data_family,
                         lb,
                         ub,
                         seed = NULL,
                         ...) {
  if (data_family == "transformed_normal") {
    link <- identity
  } else {
    link <- inv_link_lookup(data_link)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }
  dataset <- data.frame()
  iter <- 0
  while (nrow(dataset) < data_N) {
    iter <- iter + 1

    z1 <- rnorm(1, sigma_z1)
    z2 <- rnorm(1, sigma_z2)
    z3 <- rnorm(1, sigma_z3)
    x <- rnorm(1, mean = (z1_x_coef * z1 + z3_x_coef * z3), sd = sigma_x)

    mu <- do.call(
      link,
      list(
        y_intercept +
          x_y_coef * x +
          z1_y_coef * z1 +
          z2_y_coef * z2
      )
    )

    if (mu > lb & mu < ub) {
      y <- do.call(
        rng_lookup(data_family, data_link),
        list(
          1,
          mu,
          sigma_y
        )
      )

      if (y > lb & y < ub) {
        z4 <- rnorm(1, mean = (y_z4_coef * y + x_z4_coef * x), sd = sigma_z4)

        dataset <- rbind(
          dataset,
          data.frame(
            z1 = z1,
            z2 = z2,
            z3 = z3,
            z4 = z4,
            x = x,
            y = y
          )
        )
      }
    }
  }

  return(
    list(
      dataset = dataset,
      n_resample = iter - data_N
    )
  )
}


#' Title
#'
#' @param data_gen_conf
#'
#' @return
#' @export
#'
#' @examples
data_gen_conf_y_analysis <- function(data_gen_conf) {
  datagen_result <- do.call(
    basedag_data,
    data_gen_conf
  )
  dataset <- datagen_result$dataset
  n_resample <- datagen_result$n_resample
  print(paste(data_gen_conf$data_family, data_gen_conf$data_link, "data:"))
  print(paste("min:", min(dataset$y)))
  print(paste("max:", max(dataset$y)))
  print(paste("mean:", mean(dataset$y)))
  print(paste("median:", median(dataset$y)))
  print(paste("n_resample:", n_resample))
  hist(dataset$y)
}
