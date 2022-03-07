#' Title
#'
#' @param x_y_coef
#' @param ...
#' @param posterior_draws
#'
#' @return
#' @export
#'
#' @examples
p_bias <- function(posterior_draws, x_y_coef, ...) {
  return(mean(posterior_draws) - x_y_coef)
}

#' Title
#'
#' @param posterior_draws
#'
#' @return
#' @export
#'
#' @examples
p_mean <- function(posterior_draws) {
  return(mean(posterior_draws))
}

#' Title
#'
#' @param posterior_draws
#'
#' @return
#' @export
#'
#' @examples
p_sd <- function(posterior_draws) {
  return(sd(posterior_draws))
}

#' Title
#'
#' @param fit
#'
#' @return
#' @export
#'
#' @examples
p_rstar <- function(fit) {
  posterior::rstar(posterior::as_draws_array(fit))
}

#' Title
#'
#' @param posterior_draws
#' @param x_y_coef
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
rmse_s <- function(posterior_draws, x_y_coef, ...) {
  return(sqrt(mean((posterior_draws - x_y_coef)^2)))
}

#' Title
#'
#' @param posterior_draws
#' @param x_y_coef
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
mae_s <- function(posterior_draws, x_y_coef, ...) {
  return(mean(abs(posterior_draws - x_y_coef)))
}

#' Title
#'
#' @param prob
#' @param posterior_draws
#'
#' @return
#' @export
#'
#' @examples
p_quantiles <- function(posterior_draws, prob) {
  return(quantile(posterior_draws, probs = prob))
}


#' Title
#'
#' @param posterior_draws
#' @param x_y_coef
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
q_true <- function(posterior_draws, x_y_coef, ...) {
  return(length(which(posterior_draws < x_y_coef)) / length(posterior_draws))
}


#' Title
#'
#' @param fit
#'
#' @return
#' @export
#'
#' @examples
bad_pareto_ks <- function(fit) {
  psis_object <- brms:::.psis(fit, newdata = fit$data, resp = NULL)
  return(length(which(psis_object$diagnostics$pareto_k > 0.7)))
}


#' Title
#'
#' @param fit
#'
#' @return
#' @export
#'
#' @examples
sampling_time <- function(fit) {
  times <- rstan::get_elapsed_time(fit$fit)
  return(list(
    "warmup" = sum(times[, 1]),
    "sample" = sum(times[, 2])
  ))
}
