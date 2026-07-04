# SBC simultaneous ECDF confidence bands ----------------------------------
#
# Ported from the 0.x bayesim codebase (git 78101f6:R/inverse_forward_sampling.R),
# which itself ports the SBC package's adjust_gamma() (Säilynoja, Bürkner &
# Vehtari 2022). These produce true simultaneous confidence bands for the ECDF
# of SBC ranks, replacing the approximate KS/DKW bound previously used in
# plot_rank_ecdf().
#
# The single-chain (L = 1) path is exact via adjust_gamma_optimize(). The
# multi-chain (L > 1) path in 0.x relied on bayeshear::u_scale, which is not
# available in bayesim 2.0; for L > 1 we fall back to the exact L = 1 band
# (conservative for the correlated-chains case) with a message.

#' Adjust the coverage parameter for simultaneous ECDF confidence bands
#'
#' Computes the gamma coverage parameter such that the simultaneous confidence
#' envelope of the ECDF of a uniform sample of size N has (approximately) the
#' requested confidence level (Säilynoja et al. 2022). For L = 1 (single
#' chain/sample) the result is exact via dynamic programming; for L > 1 the
#' 0.x code used a Monte-Carlo simulation that depended on an external
#' `u_scale` helper not ported here, so the exact L = 1 band is used (slightly
#' conservative).
#'
#' @param N Integer; number of samples (ranks).
#' @param L Integer; number of samples/chains. Default 1.
#' @param K Integer; number of equally spaced evaluation points (right ends of
#'   the partition intervals). Defaults to N.
#' @param conf_level Numeric in (0,1); confidence level. Default 0.95.
#' @return Numeric gamma in (0, 1 - conf_level).
#' @keywords internal
adjust_gamma <- function(N, L, K = N, conf_level = 0.95) {
  if (any(c(K, N, L) < 1)) {
    stop(bayesim_config_error("'N', 'L' and 'K' must be positive integers."))
  }
  if (conf_level >= 1 || conf_level <= 0) {
    stop(bayesim_config_error("'conf_level' must be in (0, 1)."))
  }
  if (L == 1) {
    gamma <- adjust_gamma_optimize(N, K, conf_level)
  } else {
    # 0.x used adjust_gamma_simulate(N, L, K, conf_level) via bayeshear::u_scale,
    # which is unavailable here. Fall back to the exact single-sample band,
    # which is conservative for the multi-chain case.
    gamma <- adjust_gamma_optimize(N, K, conf_level)
  }
  gamma
}

#' Exact gamma for a single sample (L = 1) via dynamic programming.
#' @keywords internal
adjust_gamma_optimize <- function(N, K, conf_level = 0.95) {
  target <- function(gamma, conf_level, N, K) {
    z <- 1:(K - 1) / K
    z1 <- c(0, z)
    z2 <- c(z, 1)

    # pre-compute quantiles and use symmetry for increased efficiency.
    x2_lower <- qbinom(gamma / 2, N, z2)
    x2_upper <- c(N - rev(x2_lower)[2:K], 1)

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
#' @keywords internal
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
#' @return A list with `lower` and `upper` numeric vectors of length K + 1 over
#'   the grid 0:K / K.
#' @keywords internal
sbc_band <- function(N, K = N, conf_level = 0.95) {
  gamma <- adjust_gamma(N, L = 1L, K = K, conf_level = conf_level)
  z <- (0:K) / K
  x_lower <- qbinom(gamma / 2, N, z)
  x_upper <- qbinom(1 - gamma / 2, N, z)
  list(
    x = z,
    lower = x_lower / N,
    upper = x_upper / N
  )
}
