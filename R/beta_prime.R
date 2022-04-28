#' Probability density function for the Beta-Prime distribution (aka. inverse Beta)
#'
#' @param x Model space, defined for x >= 0, x may be scalar or vector
#' @param mu Mean parameter of pdf, mu > 0, mu must be a scalar
#' @param beta Second shape parameter of beta-prime, beta > 1 must be a scalar
#' @param log Optional argument. If true, returns log(pdf). Normally False.
#'
#' @return pdf of beta-prime, with mean parameterasation
#' @export
#'
#' @examples x <- seq(from = 0, to = 100, length.out = 1000)
#' y <- dbetaprime(x, mu = 4, beta = 2)
#' plot(x, y, type = "l", ylab = "Density", main = "dbetaprime(mu=4, beta=2)")
dbetaprime <- function(x, mu, beta, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x < 0))) {
    stop("dbetaprime is only defined for x >= 0")
  }
  if (isTRUE(beta <= 1)) {
    stop("dbetaprime is only defined for beta > 1")
  }
  if (isTRUE(mu <= 0)) {
    stop("dbetaprime is only defined for mu > 0")
  }
  if (isTRUE(length(mu) > 1 || length(beta) > 1)) {
    stop("dbetaprime only works with scalar shape-parameters")
  }


  # calculate the second argument for beta-prime, given mu
  alpha <- mu * (beta - 1)

  # alpha == 1 results in a 0 * log(x), for a possible log(0), this is 0*-Inf = NaN
  if (alpha == 1.0) {
    lpdf <- (-alpha - beta) * log1p(x) -
      log(beta(alpha, beta))
  } else  {
    # with alpha != 1, we get 0 != a <- alpha - 1. With this, we get a * -Inf,
    # which is -Inf. And for x > 0, we get a real value. So it is fixed.
    lpdf <- (alpha - 1) * log(x) +
      (-alpha - beta) * log1p(x) -
      log(beta(alpha, beta))
  }

  # return either the log or the pdf itself, given the log-value
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}


#' Title
#'
#' @param p Probabilities, for which to calculate the quantiles
#' @param mu Mean of the beta-prime, mu > 0
#' @param beta Beta Shape Argument of beta-prime, beta > 1
#'
#' @return Quantiles of the beta-prime distribution, given p, mu and beta
#' @export
#'
#' @examples x <- seq(from = 0, to = 1, length.out = 100)
#  y = bayesim::qbetaprime(x, mu = 1, beta = 2)
#  plot(x, y, type="l", ylab = "Quantile", main = "left-leaning Beta-Prime(mu=1,beta=2)"))
qbetaprime <- function(p, mu, beta) {
  # check the arguments
  if (isTRUE(any(p < 0 | p > 1))) {
    stop("p has to be in an interval of [0, 1]")
  }
  if (isTRUE(beta <= 1)) {
    stop("qbetaprime is only defined for beta > 1")
  }
  if (isTRUE(mu <= 0)) {
    stop("qbetaprime is only defined for mu > 0")
  }

  # calculate argument alpha of beta/betaprime
  alpha <- mu * (beta - 1)
  qb <- qbeta(p, alpha, beta)

  # now calculate the qbetaprime using log rules log(qbeta) - log(1 - qbeta)
  lqbp <- log(qb) - log1p(-qb)
  return(exp(lqbp))
}

#' Title
#'
#' @param n Number of beta-prime samples, to be calculated for RNG.
#' @param mu Mean of beta-prime distribution
#' @param beta Beta shape argument, of beta-prime
#'
#' @return Randomly samples distributed as beta-prime-
#' @export
#'
#' @examples y <- bayesim::rbetaprime(100, mu = 1, beta = 2)
#  hist(y, main = c(paste("Mean:", mean(y)), " for RNG of left-leaning Beta-Prime(mu=1,beta=2)"))
rbetaprime <- function(n, mu, beta) {
  # check the arguments
  if (isTRUE(beta <= 1)) {
    stop("qbetaprime is only defined for beta > 1")
  }
  if (isTRUE(mu <= 0)) {
    stop("qbetaprime is only defined for mu > 0")
  }
  return(qbetaprime(runif(n, min = 0, max = 1), mu, beta))
}

#' Title
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of betaprime given data in prep
#'
#' @examples
log_lik_betaprime <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  y <- prep$data$Y[i]
  return(dbetaprime(y, mu, beta, log = TRUE))
}


#' Title
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ...
#'
#' @return Posterior prediction of beta-prime, given data in prep
#'
#' @examples
posterior_predict_betaprime <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  return(rbetaprime(prep$ndraws, mu, beta))
}

#' Title
#'
#' @param prep BRMS data
#'
#' @return Recover the given mean of data prep
#'
#' @examples
posterior_epred_betaprime <- function(prep) {
  return(brms::get_dpar(prep, "mu"))
}


#' Title
#'
#' @param link Link function for function
#' @param link_beta Link function for beta argument
#'
#' @return BRMS beta-prime distribution family
#' @export
#'
#' @examples data <- list(a = a, y = bayesim::rbetaprime(
#'   n,
#'   exp(0.5 * rnorm(1000) + 1), 2
#' ))
#' fit1 <- brm(y ~ 1 + a,
#'   data = data, family = bayesim::betaprime(),
#'   stanvars = bayesim::betaprime()$stanvars, backend = "cmdstan"
#' )
#' plot(fit1)
betaprime <- function(link = "log", link_beta = "log1p") {
  family <- brms::custom_family(
    "betaprime",
    dpars = c("mu", "beta"),
    links = c(link, link_beta),
    lb = c(0, 1),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_betaprime,
    posterior_predict = posterior_predict_betaprime,
    posterior_epred = posterior_epred_betaprime
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real betaprime_lpdf(real y, real mu, real beta) {
        real alpha = mu * (beta - 1);
        return  ((alpha - 1) * log(y) + (-alpha - beta) * log1p(y) - log(beta(alpha, beta)));
      }

      real betaprime_rng(real mu, real beta) {
        real rb = beta_rng(mu * (beta - 1), beta);
        return (rb / (1 - rb));
      }",
    block = "functions"
  )
  return(family)
}
