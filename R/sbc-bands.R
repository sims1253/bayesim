# SBC simultaneous ECDF confidence bands ----------------------------------
#
# Ported from the 0.x bayesim codebase (git 78101f6:R/inverse_forward_sampling.R),
# which itself ports the SBC package's adjust_gamma() (Säilynoja, Bürkner &
# Vehtari 2022). The bands assume independent uniform ranks on a common support.

#' Adjust the coverage parameter for simultaneous ECDF confidence bands
#'
#' Computes the gamma coverage parameter such that the simultaneous confidence
#' envelope of the ECDF of a uniform sample of size N has (approximately) the
#' requested confidence level (Säilynoja et al. 2022), using dynamic programming
#' for independent uniform ranks on a common support.
#'
#' @param N Integer; number of samples (ranks).
#' @param K Integer; number of equally spaced evaluation points (right ends of
#'   the partition intervals). Defaults to N.
#' @param conf_level Numeric in (0,1); confidence level. Default 0.95.
#' @return Numeric gamma in (0, 1 - conf_level).
#' @noRd
adjust_gamma <- function(N, K = N, conf_level = 0.95) {
  if (
    !all(is.numeric(c(K, N))) ||
      !all(is.finite(c(K, N))) ||
      any(c(K, N) < 1) ||
      any(c(K, N) != as.integer(c(K, N)))
  ) {
    stop(bayesim_config_error("'N' and 'K' must be positive integers."))
  }
  if (
    !is.numeric(conf_level) ||
      length(conf_level) != 1L ||
      is.na(conf_level) ||
      !is.finite(conf_level) ||
      conf_level >= 1 ||
      conf_level <= 0
  ) {
    stop(bayesim_config_error("'conf_level' must be in (0, 1)."))
  }
  N <- as.integer(N)
  K <- as.integer(K)
  adjust_gamma_optimize(N, K, conf_level)
}

#' Gamma for an independent uniform sample via dynamic programming.
#' @noRd
adjust_gamma_optimize <- function(N, K, conf_level = 0.95) {
  if (K == 1L) {
    # With a single partition interval [0, 1], every empirical CDF is exactly
    # on the endpoint envelope; no dynamic-programming recursion is needed.
    return((1 - conf_level) / 2)
  }
  target <- function(gamma, conf_level, N, K) {
    z <- 1:(K - 1) / K
    z1 <- c(0, z)
    z2 <- c(z, 1)

    # pre-compute quantiles and use symmetry for increased efficiency.
    x2_lower <- qbinom(gamma / 2, N, z2)
    # seq_len(K - 1) + 1L (not 2:K) so K == 1 does not index past the vector
    # and produce NA bounds; the degenerate single-interval band stays valid.
    x2_upper <- c(N - rev(x2_lower)[seq_len(K - 1) + 1L], 1)

    # Compute the total probability of trajectories inside the confidence
    # intervals. Initialize the set and corresponding probabilities known
    # to be 0 and 1 for the starting value z1 = 0.
    x1 <- 0
    p_int <- 1
    for (i in seq_along(z1)) {
      tmp <- p_interior(
        p_int,
        x1 = x1,
        x2 = x2_lower[i]:x2_upper[i],
        z1 = z1[i],
        z2 = z2[i],
        gamma = gamma,
        N = N
      )
      x1 <- tmp$x1
      p_int <- tmp$p_int
    }
    abs(conf_level - sum(p_int))
  }
  optimize(target, c(0, 1 - conf_level), conf_level, N = N, K = K)$minimum
}

#' Interior-probability recursion helper.
#' @noRd
p_interior <- function(p_int, x1, x2, z1, z2, gamma, N) {
  z_tilde <- (z2 - z1) / (1 - z1)

  N_tilde <- rep(N - x1, each = length(x2))
  p_int <- rep(p_int, each = length(x2))
  x_diff <- outer(x2, x1, "-")
  p_x2_int <- p_int * dbinom(x_diff, N_tilde, z_tilde)

  list(p_int = rowSums(p_x2_int), x1 = x2)
}

#' Simultaneous confidence band for a uniform ECDF.
#'
#' Returns the lower/upper bounds of the simultaneous ECDF confidence envelope
#' (Säilynoja et al. 2022) at K equally spaced evaluation points, given N
#' samples and a confidence level.
#'
#' @param N Integer; number of samples (ranks).
#' @param K Integer; number of evaluation points. Defaults to N.
#' @param conf_level Numeric in (0,1); confidence level.
#' @return A list with `x` (the grid 0:K / K) and `lower` and `upper` numeric
#'   vectors of length K + 1 over that grid.
#' @noRd
sbc_band <- function(N, K = N, conf_level = 0.95) {
  gamma <- adjust_gamma(N, K = K, conf_level = conf_level)
  z <- (0:K) / K
  x_lower <- qbinom(gamma / 2, N, z)
  x_upper <- qbinom(1 - gamma / 2, N, z)
  list(
    x = z,
    lower = x_lower / N,
    upper = x_upper / N
  )
}
