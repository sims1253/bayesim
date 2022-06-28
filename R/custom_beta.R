#' Custom Beta distribibution density.
#'
#' @param x x-value, x e (0, 1)
#' @param mu Mean parameter, mean e (0, 1)
#' @param phi Precision parameter, phi > 0
#' @param log Optional argument. If TRUE, returns log(pdf). Normally False.
#'
#' @details The beta distribution has density
#' \deqn{f(y) = \frac{\Gamma(\mu\phi + (1 - \mu)\phi)}{\Gamma(\mu\phi)\Gamma((1 - \mu)\phi)} * x^{\mu\phi - 1}*(1 - x)^{(1-\mu)\phi} }
#' @details With parameterisation of the usual Beta-Distribution's shape parameters a and b as:
#' \deqn{a := \mu\phi, b := (1 - \mu)\phi}
#'
#' @return PDF of Custom Beta distribution, with mean parameterasation
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 0.99, length.out = 1000)
#' plot(x, dbeta_custom(x, mu = 0.5, phi = 1), type = "l")
dbeta_custom <- function(x, mu, phi, log = FALSE) {
  if (isTRUE(any(x <= 0 | x >= 1))) {
    stop("The value x has to be in (0, 1).")
  }
  if (isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The mean must be in (0, 1).")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("P must be above 0.")
  }
  # May be improved, by custom implementation.
  #lpdf <- dbeta(x, mu * phi, (1 - mu) * phi, log=TRUE)
  lpdf <- (log(gamma(phi)) - log(gamma(mu*phi)) - log(gamma((1 - mu)*phi))) +
   log(x) * (mu * phi - 1) + log1p(-x) * ((1 - mu) * phi - 1)
  if(log)
    return(lpdf)
  else
    return(exp(lpdf))
  # TODO: needs indepth stability testing!!!
}

#' Custom Beta distribution RNG
#'
#' @param n Number of draws.
#' @param mu Mean
#' @param phi Precision
#'
#' @return n samples beta distributed.
#' @export
#'
#' @examples hist(rbeta_custom(1000, mu=0.5, phi=1))
rbeta_custom <- function(n, mu, phi) {
  if (isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The mean must be in (0,1).")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("P must be above 0.")
  }
  return(rbeta(n, mu * phi, (1 - mu) * phi))
}

#' Log-Likelihood vignette for the Custom-Beta distribution, in Mean parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of Beta-Custom given data in prep
#'
#' @examples
log_lik_beta <- function(i, prep) {
  mu <- get_dpar(prep, "mu", i = i)
  phi <- get_dpar(prep, "phi", i = i)
  y <- prep$data$Y[i]
  return(dbeta_custom(y, mu, phi, log = TRUE))
}

#' Posterior prediction vignette for the Custom-Beta distribution, in Mean parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ...
#'
#' @return Posterior prediction of Beta-Custom, given data in prep
#'
#' @examples
posterior_predict_beta <- function(i, prep, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  phi <- get_dpar(prep, "phi", i = i)
  return(rbeta_custom(prep$ndraws, mu, phi))
}

#' Posterior expected value prediction of the custom-beta implementation.
#'
#' @param prep BRMS data
#'
#' @return Recover the given mean of data prep
#'
#' @examples
posterior_epred_beta <- function(prep) {
  mu <- get_dpar(prep, "mu")
  return(mu)
}

#' Custom-Beta BRMS-implementation in mean parametrization.
#'
#' @param link Link function for function
#' @param link_phi Link function for phi argument
#'
#' @return BRMS Beta-Custom distribution family
#'
#' @examples library(brms)
#' a <- rnorm(10000)
#' data <- list(a = a, y = bayesim::rbeta_custom(10000, bayesim::inv_logit(0.5 * a + 1), 2))
#' fit1 <- brm(y ~ 1 + a, data = data, family = bayesim::beta_custom(),
#'   stanvars = bayesim::beta_custom()$stanvars, backend = "cmdstan")
#' plot(fit1)
beta_custom <- function(link = "logit", link_phi = "log") {
  family <- brms::custom_family(
    "beta_custom",
    dpars = c("mu", "phi"),
    links = c(link, link_phi),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "real",
    log_lik = log_lik_beta,
    posterior_predict = posterior_predict_beta,
    posterior_epred = posterior_epred_beta
  )
  # TODO: optimize Stan code as well (or not! May be unused, given beta is implemented in BRMS)
  family$stanvars <- brms::stanvar(
    scode = "
      real beta_custom_lpdf(real y, real mu, real phi) {
        return(beta_lpdf(y | mu*phi, (1-mu)*phi));
      }

      real beta_custom_rng(real mu, real phi) {
        return(beta_rng(mu*phi, (1-mu)*phi));
      }",
    block = "functions"
  )
  return(family)
}
