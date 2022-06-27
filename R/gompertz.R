#' Probability density function for the Gompertz distribution, with Median parametrization.
#'
#' @param x Model space, defined for x >= 0
#' @param mu Median parameter of pdf, mu > 0
#' @param eta Second shape parameter of Gompertz, defined for eta > 1
#' @param log Optional argument. If TRUE, returns log(pdf). Normally False.
#'
#' @details PDF of Gompertz implementation, with constant b:
#' \deqn{b(\mu,\eta) := (1 / \mu) * log1p((-1 / \eta) * log(0.5))}
#' \deqn{f(x) = \eta*b*exp(\eta + bx - \eta * e^{bx})}
#'
#' @return f(x | mu, eta)
#' @export
#'
#' @examples x <- seq(from = 0, to = 10, length.out = 100)
#' plot(x, dgompertz(x, mu=2, eta=0.2), type = "l")
dgompertz <- function(x, mu, eta, log = FALSE) {
  # check arguments
  if (isTRUE(mu <= 0)) {
    stop("The Gompertz PDF is only defined for mu > 0")
  }
  if (isTRUE(eta <= 0)) {
    stop("The Gompertz PDF is only defined for eta > 0")
  }
  if (isTRUE(any(x < 0))) {
    stop("The Gompertz PDF is only defined on the positive scale")
  }

  # calculate missing argument b
  # b <- (1 / mu) * log1p((-1 / eta) * log(0.5))
  log_b <- -log(mu) + log(log1p((-1 / eta) * log(0.5)))
  b <- exp(log_b)
  # calculate log-pdf with pdf = b * eta * exp(eta + bx - eta * exp(bx))
  lpdf <- log_b + log(eta) + (eta + b * x - eta * exp(b * x))

  # now return either log or normal value
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Quantile function for the Gompertz distribution, with Median parametrization.
#'
#' @param p Quantile to be calculated
#' @param mu Median argument of Gompertz
#' @param eta Eta argument of Gompertz
#'
#' @return Inverse of CDF, calculates a value, given a probability p
#' @export
#'
#' @examples x <- seq(from = 0, to = 1, length.out = 100)
#' plot(x, qgompertz(x, mu=2, eta=0.2), type = "l")
qgompertz <- function(p, mu, eta) {
  # check arguments
  if (isTRUE(mu <= 0)) {
    stop("The Gompertz Qunatile-function is only defined for mu > 0")
  }
  if (isTRUE(eta <= 0)) {
    stop("The Gompertz Qunatile-function is only defined for eta > 0")
  }
  if (isTRUE(any(p < 0))) {
    stop("The Gompertz Qunatile-function is only defined on the positive scale")
  }

  # calculate missing argument b
  # median = (1/b) * ln[(-1/eta)*ln(1/2)+1]
  b <- (1 / mu) * log1p((-1 / eta) * log(0.5))

  # calculate the QDF as inverse of CDF
  x <- (1 / b) * log1p(-(1 / eta) * log1p(-p))

  return(x)
}

#' RNG function for the Gompertz distribution, with Median parametrization.
#'
#' @param n Number of draws
#' @param mu Median argument of Gompertz
#' @param eta Eta argument of Gompertz
#'
#' @return A Gompertz distributed RNG vector of size n
#' @export
#'
#' @examples hist(log(rgompertz(n, mu = 2, eta = 0.1)))
rgompertz <- function(n, mu, eta) {
  # check arguments
  if (isTRUE(mu <= 0)) {
    stop("The Gompertz RNG is only defined for mu > 0")
  }
  if (isTRUE(eta <= 0)) {
    stop("The Gompertz RNG is only defined for eta > 0")
  }
  return(qgompertz(runif(n, min = 0, max = 1), mu = mu, eta = eta))
}

#' Log-Likelihood vignette for the Gompertz distribution, with Median parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of gompertz given data in prep
#'
#' @examples
log_lik_gompertz <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  eta <- brms::get_dpar(prep, "eta", i = i)
  y <- prep$data$Y[i]
  return(dgompertz(y, mu, eta, log = TRUE))
}

#' Posterior-Prediction vignette for the Gompertz distribution, with Median parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ...
#'
#' @return Posterior prediction of gompertz, given data in prep
#'
#' @examples
posterior_predict_gompertz <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  eta <- brms::get_dpar(prep, "eta", i = i)
  return(rgompertz(prep$ndraws, mu, eta))
}

#' Expectation-Predict vignette for the Gompertz distribution, with Median parametrization.
#' Not defined for the Gompertz family.
#'
#' @param prep BRMS data
#'
#' @return Recover the given mean of data prep
#'
#' @examples
posterior_epred_gompertz <- function(prep) {
  stop("posterior_epred is not defined for the gompertz family")
}


#' Custom Gompertz BRMS-implementation in median parametrization.
#'
#' @param link Link function for function
#' @param link_eta Link function for eta argument
#'
#' @return BRMS gompertz distribution family
#' @export
#'
#' @examples a <- rnorm(10000)
#' data <- list(a = a, y = bayesim::rgompertz(10000, exp(0.5 * a + 1), 0.2))
#' fit1 <- brm(y ~ 1 + a, data = data, family = bayesim::gompertz(),
#'   stanvars = bayesim::gompertz()$stanvars, backend = "cmdstan")
#' plot(fit1)
gompertz <- function(link = "log", link_eta = "log") {
  family <- brms::custom_family(
    "gompertz",
    dpars = c("mu", "eta"),
    links = c(link, link_eta),
    lb = c(0, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_gompertz,
    posterior_predict = posterior_predict_gompertz,
    posterior_epred = posterior_epred_gompertz
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real gompertz_lpdf(real y, real mu, real eta) {
        real log_b = -log(mu) + log(log1p((-1 / eta) * log(0.5)));
        real b = exp(log_b);
        real lpdf = log_b + log(eta) + (eta + b * y - eta * exp(b * y));
        return(lpdf);
      }
      real gompertz_rng(real mu, real eta) {
        return((1 / ((1 / mu) * log1p(-(1 / eta) * log(0.5)))) *
                log1p((-1 / eta) *
                log(uniform_rng(0, 1))));
      }",
    block = "functions"
  )
  return(family)
}
