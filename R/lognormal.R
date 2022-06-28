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
#' @examples x <- seq(from = 0.01, to = 10, length.out = 1000)
#' plot(x, dlognormal(x, mu = 2, sigma = 2), type = "l")
dlognormal <- function(x, mu, sigma, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("lognormal is only defined for x > 0")
  }
  if (isTRUE(sigma <= 0)) {
    stop("lognormal is only defined for sigma > 0")
  }
  logpdf <-
    -(log(x) + log(sigma) + 0.5 * (log(2) + log(pi))) +
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
#' @examples hist(log(rlognormal(100, 0.5, 2)))
rlognormal <- function(n, mu, sigma) {
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
#' @return log_likelihood of the Lognormal distribution, given some BRMS data.
log_lik_lognormal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dlognormal(y, mu, sigma, log = TRUE))
}

#' Posterior-predict vignette for the Lognormal distribution, with Median parametrization.
#'
#' @param i Indices
#' @param prep BRMS data
#' @param ...
#'
#' @return The posterior prediction of the Lognormal distribution, given some BRMS data.
posterior_predict_lognormal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rlognormal(prep$ndraws, mu, sigma))
}

#' Posterior expected value prediction vignette for Lognormal distribution.
#'
#' @param prep BRMS data
#'
#' @return Mean of Posterior
posterior_epred_lognormal <- function(prep) {
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
#' @examples library(brms)
#' a <- rnorm(1000)
#' data <- list(a = a, y = rlognormal(n, exp(0.5 * a + 1), 2))
#' fit1 <- brm(y ~ 1 + a, data = data, family = lognormal(),
#'   stanvars = lognormal()$stanvars, backend = "cmdstan")
#' plot(fit1)
lognormal <- function(link = "identity", link_sigma = "log") {
  stopifnot(link == "identity")
  family <- brms::custom_family(
    "lognormal",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_lognormal,
    posterior_predict = posterior_predict_lognormal,
    posterior_epred = posterior_epred_lognormal
  )
  family$stanvars <- stanvars <- brms::stanvar(
    scode = "
      real lognormal_lpdf(real y, real mu, real sigma) {
        return -(log(y) + log(sigma) + 0.5 * (log(2) + log(pi()))) +
                (-(log(y) - mu)^2 / (2 * sigma^2));
      }

      real lognormal_rng(real mu, real sigma) {
        return exp(normal_rng(mu, sigma));
      }",
    block = "functions"
  )
  return(family)
}
