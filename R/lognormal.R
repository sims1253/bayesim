#' Lognormal density distribution in median parametrization.
#'
#' @param x Value space of the distribution, x > 0
#' @param mu Median parameter, mu is already log-transformed, mu unbound
#' @param sigma Sigma shape parameter, sigma >= 0
#' @param log Bool argument, if true, returns the logarithmic density
#'
#' @return Normal distribution density with logit link function
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 10, length.out = 100)
#' plot(x, dlognormal_custom(x, mu = 1, sigma = 0.5), type = "l")
dlognormal_custom <- function(x, mu, sigma, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("lognormal is only defined for x > 0")
  }
  if (isTRUE(sigma <= 0)) {
    stop("lognormal is only defined for sigma > 0")
  }
  logpdf <-
    -(log(sigma) + 0.5 * (log(2 * pi))) +
    -log(x) +
    (-(log(x) - mu)^2 / (2 * sigma^2))
  if (log) {
    return(logpdf)
  } else {
    return(exp(logpdf))
  }
}

#' Lognormal RNG-function in median parametrization.
#'
#' @param n Number of draws
#' @param mu Median paramameter, mu unbound, mu already log transformed
#' @param sigma Sigma shape parameter, sigma > 0
#'
#' @returns n Lognormally ditributed samples
#'
#' @export
#'
#' @examples hist(rlognormal_custom(100, 1, 0.5))
rlognormal_custom <- function(n, mu, sigma) {
  # check the arguments
  if (isTRUE(sigma <= 0)) {
    stop("lognormal is only defined for sigma > 0")
  }
  return(
    exp(rnorm(n, mu, sigma))
  )
}

#' Log-Likelihood vignette for the Lognormal distribution, in Median parametrization.
#'
#' @param i Indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood for BRMS
#'
log_lik_lognormal_custom <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dlognormal_custom(y, mu, sigma, log = TRUE))
}

#' Posterior-predict vignette for the Lognormal distribution, with Median parametrization.
#'
#' @param i Indices
#' @param prep BRMS data
#' @param ...
#'
#' @return Posterior prediction for BRMS
#'
posterior_predict_lognormal_custom <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rlognormal_custom(prep$ndraws, mu, sigma))
}

#' Posterior expected value prediction vignette for Lognormal distribution.
#'
#' @param prep BRMS data
#'
#' @return Expectation value prediction for BRMS
#'
posterior_epred_lognormal_custom <- function(prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(exp(mu + sigma^2 / 2))
}

#' Custom BRMS family Log-Normal in median parametrization.
#'
#' @param link Link function argument (as string) for Median argument. Left as identity!
#' @param link_sigma Link function argument (as string) for Shape argument
#'
#' @return Lognormal BRMS model-object
#' @export
#'
#' @examples # Running the example might take a while and may make RStudio unresponsive.
#' # Just relax and grab a cup of coffe or tea in the meantime.
#' library(bayesim)
#' library(BBmisc)
#' library(brms)
#' a <- rnorm(1000)
#' data <- list(a = a, y = bayesim::rlognormal_custom(1000, 0.5 * a + 1, 2))
#' # BBmisc::surpressAll necassary, the RStudio Roxygen help would be filled with slash symbols...
#' # For an example without surpress, checkout the Bayesim Betaprime Example script
#' BBmisc::suppressAll({  fit1 <- brms::brm(y ~ 1 + a, data = data, family = bayesim::lognormal_custom(),
#'   stanvars = bayesim::lognormal_custom()$stanvars, backend = "cmdstanr", cores = 4)  })
#' plot(fit1)
lognormal_custom <- function(link = "identity", link_sigma = "log") {
  stopifnot(link == "identity")
  family <- brms::custom_family(
    "lognormal_custom",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_lognormal_custom,
    posterior_predict = posterior_predict_lognormal_custom,
    posterior_epred = posterior_epred_lognormal_custom
  )
  family$stanvars <- stanvars <- brms::stanvar(
    scode = "
      real lognormal_custom_lpdf(real y, real mu, real sigma) {
        return -(log(y) + log(sigma) + 0.5 * (log(2 * pi()))) +
               -(log(y) - mu)^2 / (2 * sigma^2);
      }

      real lognormal_custom_rng(real mu, real sigma) {
        return exp(normal_rng(mu, sigma));
      }",
    block = "functions"
  )
  return(family)
}
