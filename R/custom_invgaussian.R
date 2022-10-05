#' Probability density function for the Inverse-Gaussion distribution
#'
#' @param x Value space, x > 0.
#' @param mu Mean parameter of the density, mu > 0.
#' @param shape shape > 0
#'
#' @details The PDF is defined with shape parameter alpha as
#' \deqn{f(y) = \sqrt{\frac{\alpha}{2 \pi y^3}} exp(\frac{-\alpha(y-\mu)^2}{2 \mu^2 y}) }
#' \url{https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution}
#'
#' @return f(x | mu, shape)
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 10, length.out = 1000)
#' plot(x, dinversegaussian_custom(x, mu = 2, shape = 2), type = "l")
dinversegaussian_custom <- function(x, mu, shape, log = FALSE) {
  if (isTRUE(any(x <= 0))) {
    stop("Inverse Gaussian density is only defined for x > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("Inverse Gaussian density is only defined for mu > 0")
  }
  if (isTRUE(shape <= 0)) {
    stop("Inverse Gaussian density is only defined for shape > 0")
  }
  lpdf <- 0.5 * (log(shape) - (log(2) + log(pi) + 3 * log(x))) -
    shape * (x - mu)^2 / (2 * mu^2 * x)
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Custom Inverse-Gaussian RNG function in mean parametrisation.
#'
#' @param n Number of samples, scalar natural number.
#' @param mu Mean parameter, mu > 0.
#' @param shape shape > 0.
#'
#' @return n Inverse-Gaussian distributed samples.
#' @export
#'
#' @examples hist(log(rinversegaussian_custom(10000, mu = 2, shape = 1)))
rinversegaussian_custom <- function(n, mu, shape) {
  if (isTRUE(mu <= 0)) {
    stop("Inverse Gaussian density is only defined for mu > 0")
  }
  if (isTRUE(shape <= 0)) {
    stop("Inverse Gaussian density is only defined for shape > 0")
  }
  brms::rinv_gaussian(n, mu, shape)
}


#' Log-Likelihood vignette for the Inverse-Gaussian distribution, in mean parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of inversegaussian given data in prep
log_lik_inversegaussian_custom <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  shape <- brms::get_dpar(prep, "shape", i = i)
  y <- prep$data$Y[i]
  return(dinversegaussian_custom(y, mu, shape, log = TRUE))
}


#' Posterior-Predict vignette for the Inverse-Gaussian distribution, in mean parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ...
#'
#' @return Posterior prediction of inversegaussian, given data in prep
posterior_predict_inversegaussian_custom <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  shape <- brms::get_dpar(prep, "shape", i = i)
  return(rgompertz(prep$ndraws, mu, shape))
}


#' Posterior-Expectation-Predict vignette for the Inverse-Gaussian distribution, in mean parametrization.
#'
#' @param prep BRMS data
#'
#' @return Recover the given mean of data prep
posterior_epred_inversegaussian_custom <- function(prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  return(mu)
}



#' Custom Inverse-Gaussian BRMS-implementation in mean parametrization.
#'
#' @param link Link function for function
#' @param link_shape Link function for the shape argument
#'
#' @return BRMS inversegaussian distribution family
#' @export
#'
#' @examples # Running the example might take a while and may make RStudio unresponsive.
#' # Just relax and grab a cup of coffe or tea in the meantime.
#' a <- rnorm(1000)
#' data <- list(a = a, y = bayesim::rinversegaussian_custom(1000, exp(0.5 * a + 1), 2))
#' # BBmisc::surpressAll necassary, the RStudio Roxygen help would be filled with slash symbols...
#' # For an example without surpress, checkout the Bayesim Betaprime Example script
#' BBmisc::suppressAll({
#'   fit1 <- brms::brm(y ~ 1 + a,
#'     data = data, family = bayesim::inversegaussian_custom(),
#'     stanvars = bayesim::inversegaussian_custom()$stanvars, backend = "cmdstanr", cores = 4
#'   )
#' })
#' plot(fit1)
inversegaussian_custom <- function(link = "log", link_shape = "log") {
  family <- brms::custom_family(
    "inversegaussian_custom",
    dpars = c("mu", "shape"),
    links = c(link, link_shape),
    lb = c(0, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_inversegaussian_custom,
    posterior_predict = posterior_predict_inversegaussian_custom,
    posterior_epred = posterior_epred_inversegaussian_custom
  )
  family$stanvars <- brms::stanvar(
    scode = "
        real inversegaussian_custom_lpdf(real y, real mu, real shape) {
          return (
            0.5 *(log(shape) - (log(2) + log(pi()) + 3 * log(y))) -
            shape * (y - mu)^2/(2 * mu^2 * y)
          );
        }",
    block = "functions"
  )
  return(family)
}
