#' Title
#'
#' @param x
#' @param mu
#' @param phi
#' @param log
#'
#' @return
#' @export
#'
#' @examples
dbeta_custom <- function(x, mu, phi, log = FALSE) {
  if (isTRUE(any(x <= 0 || x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(mu <= 0 || mu >= 1))) {
    stop("The mean must be in (0,1).")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("phi must be above 0.")
  }
  res <- (mu * phi -1) * log(x) +
    ((1-mu )* phi -1) * log1p(-x) -
    lbeta(mu * phi, (1 - mu) * phi)
  if(log) {
    return(res)
  } else {
    return(exp(res))
  }
}

#' Title
#'
#' @param n
#' @param mu
#' @param phi
#'
#' @return
#' @export
#'
#' @examples
rbeta_custom <- function(n, mu, phi) {
  if (isTRUE(any(mu <= 0 || mu >= 1))) {
    stop("The median must be in (0,1).")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("P must be above 0.")
  }
  return(rbeta(n, mu * phi, (1 - mu) * phi))
}


#' Title
#'
#' @param i
#' @param prep
#'
#' @return
#'
#' @examples
log_lik_beta_custom <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  y <- prep$data$Y[i]
  return(dbeta_custom(y, mu, phi, log = TRUE))
}

#' Title
#'
#' @param i
#' @param prep
#' @param ...
#'
#' @return
#'
#' @examples
posterior_predict_beta_custom <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  p <- brms::get_dpar(prep, "phi", i = i)
  return(rbeta_custom(prep$ndraws, mu, phi))
}


#' Title
#'
#' @param prep
#'
#' @return
#'
#' @examples
posterior_epred_beta_custom <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  return(mu)
}

#' Title
#'
#' @param link
#' @param link_phi
#'
#' @return
#' @export
#'
#' @examples
beta_custom <- function(link = "logit", link_phi = "log") {
  family <- brms::custom_family(
    "beta_custom",
    dpars = c("mu", "phi"),
    links = c(link, link_phi),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "real",
    log_lik = log_lik_beta_custom,
    posterior_predict = posterior_predict_beta_custom,
    posterior_epred = posterior_epred_beta_custom
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real beta_custom_lpdf(real y, real mu, real phi) {
        return (mu * phi -1) * log(y) +
               ((1-mu )* phi -1) * log1m(y) -
               lbeta(mu * phi, (1 - mu) * phi);
    }

      real beta_custom_rng(real mu, real phi) {
         return beta_custom_rng(mu * phi, (1 - mu) * phi)
      }",
    block = "functions"
  )
  return(family)
}
