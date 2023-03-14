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
#' @param data_family
#' @param lb
#' @param ub
#' @param resample
#' @param seed
#' @param testing_data
#' @param ...
#' @param noise_sd
#' @param noisy_z1_x_coef
#' @param noisy_z3_x_coef
#' @param noisy_x_y_coef
#' @param noisy_z1_y_coef
#' @param noisy_z2_y_coef
#' @param noisy_y_z4_coef
#' @param noisy_x_z4_coef
#' @param noisy_y_intercept
#' @param noisy_sigma_z1
#' @param noisy_sigma_z2
#' @param noisy_sigma_z3
#' @param noisy_sigma_z4
#' @param noisy_sigma_x
#' @param noisy_sigma_y
#'
#' @return
#' @export
#'
#' @examples
basedag_data_noisy <- function(data_N,
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
                               noise_sd = 1,
                               noisy_z1_x_coef = NULL,
                               noisy_z3_x_coef = NULL,
                               noisy_x_y_coef = NULL,
                               noisy_z1_y_coef = NULL,
                               noisy_z2_y_coef = NULL,
                               noisy_y_z4_coef = NULL,
                               noisy_x_z4_coef = NULL,
                               noisy_y_intercept = NULL,
                               noisy_sigma_z1 = NULL,
                               noisy_sigma_z2 = NULL,
                               noisy_sigma_z3 = NULL,
                               noisy_sigma_z4 = NULL,
                               noisy_sigma_x = NULL,
                               noisy_sigma_y = NULL,
                               ...) {
  noisy_list <- list(
    "noisy_z1_x_coef",
    "noisy_z3_x_coef",
    "noisy_x_y_coef",
    "noisy_z1_y_coef",
    "noisy_z2_y_coef",
    "noisy_y_z4_coef",
    "noisy_x_z4_coef",
    "noisy_y_intercept",
    "noisy_sigma_z1",
    "noisy_sigma_z2",
    "noisy_sigma_z3",
    "noisy_sigma_z4",
    "noisy_sigma_x"
  )
  for (argument in noisy_list) {
    if (is.null(get(argument))) {
      assign(
        argument,
        rnorm(1,
          mean = get(substring(argument, first = 7)),
          sd = noise_sd * get(substring(argument, first = 7))
        )
      )
    }
  }
  if (is.null(noisy_sigma_y)) {
    bounds <- aux_limits_lookup(data_family)
    while (is.null(noisy_sigma_y)) {
      sample_list <- rnorm(100, mean = sigma_y, sd = noise_sd * sigma_y)
      sample_list <- subset(
        sample_list,
        sample_list > bounds$lb & sample_list < bounds$ub
      )
      if (length(sample_list) > 0) {
        noisy_sigma_y <- sample_list[[1]]
      }
    }
  }
  arguments <- as.list(c(as.list(environment()), list(...)))
  arguments$seed <- NULL
  arguments$noisy_list <- NULL
  arguments$argument <- NULL
  arguments$bounds <- NULL
  arguments$sample_list <- NULL

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
    z1 <- rnorm(n = resample * data_gen_size, mean = 0, sd = noisy_sigma_z1)
    z2 <- rnorm(n = resample * data_gen_size, mean = 0, sd = noisy_sigma_z2)
    z3 <- rnorm(n = resample * data_gen_size, mean = 0, sd = noisy_sigma_z3)
    x <- rnorm(n = resample * data_gen_size, mean = (noisy_z1_x_coef * z1 + noisy_z3_x_coef * z3), sd = noisy_sigma_x)

    mu <- do.call(
      link_lookup(data_link, inv = TRUE),
      list(
        y_intercept +
          noisy_x_y_coef * x +
          noisy_z1_y_coef * z1 +
          noisy_z2_y_coef * z2
      )
    )

    y <- do.call(
      rng_lookup(data_family),
      list(
        length(mu),
        mu,
        noisy_sigma_y
      )
    )

    sanitized_index <- which(y < ub & y > lb)
    z1 <- z1[sanitized_index]
    z2 <- z2[sanitized_index]
    z3 <- z3[sanitized_index]
    x <- x[sanitized_index]
    y <- y[sanitized_index]

    z4 <- rnorm(n = length(y), mean = (noisy_y_z4_coef * y + noisy_x_z4_coef * x), sd = noisy_sigma_z4)

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

  data_gen_output <- list(
    sampling_loops = iter,
    bad_samples = bad_samples
  )
  data_gen_output <- c(data_gen_output, arguments)

  if (testing_data) {
    return(
      list(
        dataset = dataset[1:data_N, ],
        testing_data = dataset[(data_N + 1):data_gen_size, ],
        data_gen_output = data_gen_output
      )
    )
  } else {
    return(
      list(
        dataset = dataset,
        testing_data = NULL,
        data_gen_output = data_gen_output
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
