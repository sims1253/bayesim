#' Simulate a new dataset using forward sampling.
#'
#' @param fit An object of class brmsfit or a list of brmsfit objects.
#' @param i The index of a single posterior draw to simulate a dataset for.
#'  The index is passed to [posterior_predict()]'s "draw_ids"
#'  argument.
#' @param newdata A dataframe that is passed to posterior predict.
#' @param ... Potential additional arguments.
#'
#' @return A data.frame containing n observations for each variable in the fit.
#' @export forward_sampling forward_sampling.list
#'
#'
forward_sampling <- function(fit, i, n, ...) {
  UseMethod("forward_sampling")
}

#'
forward_sampling.list <- function(x, i, n, ...) {
  if (any(lapply(x, bayeshear::post_warmup_samples) < i) | i <= 0) {
    stop("You tried to use a non-existent posterior sample.")
  }

  resp_list <- list()
  var_list <- list()

  for (index in seq_along(x)) {
    fit <- x[[index]]
    if (is(fit$formula, "brmsformula")) {
      resp_list[[fit$formula$resp]] <- index
      var_list[[fit$formula$resp]] <- all.vars(fit$formula$formula)[-1]
    } else if (is(fit$formula, "mvbrmsformula")) {
      for (formula in fit$formula$forms) {
        resp_list[[formula$resp]] <- index
        var_list[[formula$resp]] <- all.vars(formula$formula)[-1]
      }
    } else {
      stop(
        "Unsupported model type detected!
           Please use brmsfit or stanfit objects."
      )
    }
  }
  df <- data.frame(matrix(
    0, # Create empty data frame
    nrow = n,
    ncol = length(resp_list)
  ))
  colnames(df) <- names(resp_list)

  missing_responses <- names(resp_list)
  while (length(missing_responses) > 0) {
    success <- FALSE
    for (response in missing_responses) {
      if (length(intersect(missing_responses, var_list[[response]])) == 0) {
        df[response] <- as.vector(
          brms::posterior_predict(
            x[[resp_list[[response]]]],
            newdata = df[var_list[[response]]],
            resp = response,
            draw_ids = i,
            ...
          )
        )
        missing_responses <- missing_responses[!(missing_responses == response)]
        success <- TRUE
      }
    }
    if (!success) {
      stop("Unsolvable variable dependencies detected.")
    }
  }
  return(df)
}

#' Full forward sampling of a the response of brms fit, including multivariate models.
#'
#' @source This function is taken from the SBC package
#'   (https://github.com/hyunjimoon/SBC) and only here to ensure stability, as
#'   it is not exported in SBC.
#'
#' @param fit An object of class `brmsfit`
#' @param newdata An optional data.frame for which to evaluate predictions. If NULL (default), the original data of the model is used.
#' @param draws An integer vector specifying the posterior draws to be used. If NULL (the default), all draws are used.
#' @param validate_all if TRUE, validation of input data will be done in all iterations, otherwise only once
#'
#' @return A list of data.frames containing the draws.
#'
#' @keywords internal
#' @examples # Pending
#' @export
brms_full_ppred <- function(
  fit,
  newdata = NULL,
  draws = NULL,
  validate_all = FALSE
) {
  # 1. determine term hierarchy
  resp <- brms_response_sequence(fit)
  # 2.1. initialize dataframe using the original fit's data
  if (is.null(newdata)) {
    newdata <- fit$data
  }
  n <- nrow(newdata)
  # 2.3. if no draws set, range from 1 to all iters (check draws < iters)
  if (is.null(draws)) {
    draws <- seq_len(sum(fit$fit@sim$n_save))
  }
  # 2.4. create list to hold data
  pp_data <- list()

  if (!validate_all) {
    # Validate once
    newdata <- brms::validate_newdata(newdata, fit, allow_new_levels = TRUE)
  }

  for (i in draws) {
    pp_data[[i]] <- newdata
    for (vars in resp) {
      pp_data[[i]][, vars] <- array(
        brms::posterior_predict(
          fit,
          newdata = pp_data[[i]],
          resp = vars,
          draw_ids = i,
          skip_validate = !validate_all,
          allow_new_levels = TRUE
        ),
        dim = c(1, n, length(vars))
      )[1, , ]
    }
  }
  pp_data
}

nodes_by_depth <- function(adj_matrix) {
  depth_list <- list()
  var_names <- rownames(adj_matrix)
  while (nrow(adj_matrix)) {
    pos <- which(apply(adj_matrix, 1, sum) == 0)
    depth_list <- c(depth_list, list(var_names[pos]))

    var_names <- var_names[-pos]
    adj_matrix <- adj_matrix[-pos, -pos, drop = FALSE]
  }
  depth_list
}

#' Determine the response sequence of brms model
#' @source This function is taken from the SBC package
#'   (https://github.com/hyunjimoon/SBC) and only here to ensure stability, as
#'   it is not exported in SBC.
#' @keywords internal
brms_response_sequence <- function(x) {
  UseMethod("brms_response_sequence")
}

#' @source This function is taken from the SBC package
#'   (https://github.com/hyunjimoon/SBC) and only here to ensure stability, as
#'   it is not exported in SBC.
#' @method brms_response_sequence brmsfit
#' @export
brms_response_sequence.brmsfit <- function(x) {
  brms_response_sequence(x$formula)
}

#' @source This function is taken from the SBC package
#'   (https://github.com/hyunjimoon/SBC) and only here to ensure stability, as
#'   it is not exported in SBC.
#' @method brms_response_sequence bform
#' @export
brms_response_sequence.bform <- function(x) {
  term_list <- brms_response_sequence(brms::brmsterms(x))
  resp_vars <- names(term_list)

  adjacency <- t(sapply(term_list, \(x) is.element(resp_vars, x)))
  attr(adjacency, "dimnames") <- list(resp_vars, resp_vars)
  nodes_by_depth(adjacency)
}

#' @source This function is taken from the SBC package
#'   (https://github.com/hyunjimoon/SBC) and only here to ensure stability, as
#'   it is not exported in SBC.
#' @method brms_response_sequence mvbrmsterms
#' @export
brms_response_sequence.mvbrmsterms <- function(x) {
  names(x$terms) <- NULL
  sapply(x$terms, brms_response_sequence)
}

#' @source This function is taken from the SBC package
#'   (https://github.com/hyunjimoon/SBC) and only here to ensure stability, as
#'   it is not exported in SBC.
#' @method brms_response_sequence brmsterms
#' @export
brms_response_sequence.brmsterms <- function(x) {
  vars <- list(unique(unlist(lapply(x$dpars, brms_response_sequence))))
  names(vars) <- all.vars(x$respform)
  vars
}

#' @source This function is taken from the SBC package
#'   (https://github.com/hyunjimoon/SBC) and only here to ensure stability, as
#'   it is not exported in SBC.
#' @method brms_response_sequence btl
#' @export
brms_response_sequence.btl <- function(x) {
  c("1", all.vars(x$formula))
}


#' @export
#' @param N - length of samples (chains).
#' @param L - number of samples (chains).
#' @param K - number of equally spaced evaluation points, i.e. the right ends of the
#' partition intervals.
adjust_gamma <- function(N, L, K = N, conf_level = 0.95) {
  if (any(c(K, N, L) < 1)) {
    abort("Parameters 'N', 'L' and 'K' must be positive integers.")
  }
  if (conf_level >= 1 || conf_level <= 0) {
    abort("Value of 'conf_level' must be in (0,1).")
  }
  if (L == 1) {
    gamma <- adjust_gamma_optimize(N, K, conf_level)
  } else {
    gamma <- adjust_gamma_simulate(N, L, K, conf_level)
  }
  gamma
}

# Adjust coverage parameter to find silmultaneous confidence intervals for the
# ECDF of a sample from the uniform distribution.
# N - length of samples
# K - number of equally spaced evaluation points, i.e. the right ends of the
# partition intervals.
adjust_gamma_optimize <- function(N, K, conf_level = 0.95) {
  target <- function(gamma, conf_level, N, K) {
    z <- 1:(K - 1) / K
    z1 <- c(0, z)
    z2 <- c(z, 1)

    # pre-compute quantiles and use symmetry for increased efficiency.
    x2_lower <- qbinom(gamma / 2, N, z2)
    x2_upper <- c(N - rev(x2_lower)[2:K], 1)

    # Compute the total probability of trajectories inside the confidence
    # intervals. Initialize the set and corresponding probasbilities known
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

# Adjust coverage parameter to find silmultaneous confidence intervals for the
# ECDFs of multiple samples (chains) from the uniform distribution.
# N - length of samples (chains).
# L - number of samples (chains).
# K - number of equally spaced evaluation points, i.e. the right ends of the
# partition intervals.
# M - number of simulations used to determine the 'conf_level' middle quantile.
adjust_gamma_simulate <- function(N, L, K, conf_level = 0.95, M = 5000) {
  gamma <- numeric(M)
  z <- (1:(K - 1)) / K
  if (L > 1) {
    n <- N * (L - 1)
    k <- floor(z * N * L)
    for (m in seq_len(M)) {
      u <- u_scale(replicate(L, runif(N)))
      scaled_ecdfs <- apply(outer(u, z, "<="), c(2, 3), sum)
      gamma[m] <- 2 *
        min(
          apply(
            scaled_ecdfs,
            1,
            phyper,
            m = N,
            n = n,
            k = k
          ),
          apply(
            scaled_ecdfs - 1,
            1,
            phyper,
            m = N,
            n = n,
            k = k,
            lower.tail = FALSE
          )
        )
    }
  } else {
    for (m in seq_len(M)) {
      u <- runif(N)
      scaled_ecdf <- colSums(outer(u, z, "<="))
      gamma[m] <- 2 *
        min(
          pbinom(scaled_ecdf, N, z),
          pbinom(scaled_ecdf - 1, N, z, lower.tail = FALSE)
        )
    }
  }
  alpha_quantile(gamma, 1 - conf_level)
}

p_interior <- function(p_int, x1, x2, z1, z2, gamma, N) {
  z_tilde <- (z2 - z1) / (1 - z1)

  N_tilde <- rep(N - x1, each = length(x2))
  p_int <- rep(p_int, each = length(x2))
  x_diff <- outer(x2, x1, "-")
  p_x2_int <- p_int * dbinom(x_diff, N_tilde, z_tilde)

  list(p_int = rowSums(p_x2_int), x1 = x2)
}

# 100 * `alpha` percent of the trials are allowed to be rejected.
# In case of ties, return the largest value dominating at most
# 100 * (alpha + tol) percent of the values.
alpha_quantile <- function(gamma, alpha, tol = 0.001) {
  a <- unname(quantile(gamma, probs = alpha))
  a_tol <- unname(quantile(gamma, probs = alpha + tol))
  if (a == a_tol) {
    if (min(gamma) < a) {
      # take the largest value that doesn't exceed the tolerance.
      a <- max(gamma[gamma < a])
    }
  }
  a
}

#'  Quantifies deviation from uniformity by the likelihood of observing the
#'  most extreme point on the empirical CDF of the given rank distribution
#'  according to [1] (equation 7).
#'
#'  [1] Modrák, Martin, Angie H. Moon, Shinyoung Kim, Paul Bürkner, Niko Huurre,
#'  Kateřina Faltejsková, Andrew Gelman, and Aki Vehtari.
#'  “Simulation-Based Calibration Checking for Bayesian Computation:
#'  The Choice of Test Quantities Shapes Sensitivity.” arXiv, June 15, 2023.
#'  https://doi.org/10.48550/arXiv.2211.02383.
#'
#' @param ranks Rank distribution
#' @param post_warmup_draws Number of posterior draws that were used to
#'  calculate the rank distribution.
#' @param log True of the result should be on the log scale.
#'
#' @return Measure quantifying deviation from uniformity. This value can
#'  be compared to the distribution of gamma expected under uniformity
#'  calculated by validation.gamma_null_distribution.
#' @export
#'
#' @examples
gamma_discrepancy <- function(ranks, post_warmup_draws, log = FALSE) {
  if (any(is.na(ranks))) {
    return(NA)
  }
  # observed count of ranks smaller than i
  R_i <- sapply(1:(post_warmup_draws + 1), function(i) sum(ranks < i))
  # expected proportion of observed ranks smaller than i
  z_i <- sapply(
    1:(post_warmup_draws + 1),
    function(i) i / (post_warmup_draws + 1)
  )

  x1 <- pbinom(q = R_i, size = length(ranks), prob = z_i)
  x2 <- 1 - pbinom(q = R_i - 1, size = length(ranks), prob = z_i)

  if (log) {
    return(log(2 * min(x1, x2, na.rm = TRUE)))
  } else {
    return(2 * min(x1, x2, na.rm = TRUE))
  }
}
