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
#' @param resample
#' @param testing_data
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
                         resample = 1.3,
                         seed = NULL,
                         testing_data = TRUE,
                         ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  dataset <- data.frame()
  iter <- 0
  if (testing_data) {
    data_gen_size <- data_N * 2
  } else {
    data_gen_size <- data_N
  }

  while (nrow(dataset) < data_gen_size) {
    z1 <- rnorm(n = resample * data_gen_size, mean = 0, sd = sigma_z1)
    z2 <- rnorm(n = resample * data_gen_size, mean = 0, sd = sigma_z2)
    z3 <- rnorm(n = resample * data_gen_size, mean = 0, sd = sigma_z3)
    x <- rnorm(n = resample * data_gen_size, mean = (z1_x_coef * z1 + z3_x_coef * z3), sd = sigma_x)

    mu <- do.call(
      inv_link_lookup(data_link),
      list(
        y_intercept +
          x_y_coef * x +
          z1_y_coef * z1 +
          z2_y_coef * z2
      )
    )

    y <- do.call(
      rng_lookup(data_family),
      list(
        length(mu),
        mu,
        sigma_y
      )
    )

    sanitized_index <- which(y < ub & y > lb)
    z1 <- z1[sanitized_index]
    z2 <- z2[sanitized_index]
    z3 <- z3[sanitized_index]
    x <- x[sanitized_index]
    y <- y[sanitized_index]

    z4 <- rnorm(n = length(y), mean = (y_z4_coef * y + x_z4_coef * x), sd = sigma_z4)

    dataset <- rbind(
      dataset,
      data.frame(
        z1 = z1[1:min((data_gen_size - nrow(dataset)), length(y))],
        z2 = z2[1:min((data_gen_size - nrow(dataset)), length(y))],
        z3 = z3[1:min((data_gen_size - nrow(dataset)), length(y))],
        z4 = z4[1:min((data_gen_size - nrow(dataset)), length(y))],
        x = x[1:min((data_gen_size - nrow(dataset)), length(y))],
        y = y[1:min((data_gen_size - nrow(dataset)), length(y))]
      )
    )
    iter <- iter + 1
    bad_samples <- (data_gen_size * resample) - length(y)
  }

  if (testing_data) {
    return(
      list(
        dataset = dataset[1:data_N, ],
        testing_data = dataset[(data_N + 1):data_gen_size, ],
        sampling_loops = iter,
        bad_samples = bad_samples
      )
    )
  } else {
    return(
      list(
        dataset = dataset,
        testing_data = NULL,
        sampling_loops = iter,
        bad_samples = bad_samples
      )
    )
  }
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
  print(paste("y_intercept:", data_gen_conf$y_intercept, "y_sigma:", data_gen_conf$sigma_y))
  print(paste("min:", min(dataset$y)))
  print(paste("max:", max(dataset$y)))
  print(paste("mean:", mean(dataset$y)))
  print(paste("median:", median(dataset$y)))
  print(paste("n_resample:", n_resample))
  hist(dataset$y, main = paste(data_gen_conf$data_family, data_gen_conf$data_link, data_gen_conf$shape, "x_y_coef: ", data_gen_conf$x_y_coef), breaks = 20)
}
