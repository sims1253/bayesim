# Built-in study parts --------------------------------------------------------------
#
# A generator that draws the truth from a prior, an analytic method, four
# measurements of parameter recovery and calibration, and one paired
# comparison. Together they make a complete study without a sampler.

#' Generate data from parameters drawn from a prior
#'
#' `prior_draws()` makes a [study()] generator. For each dataset it draws the
#' true parameters with `prior` and then the data with `simulate`. When the
#' methods fit the same prior, the posterior ranks of the truth are uniform;
#' [measure_sbc_rank()] checks this.
#'
#' @param prior `function(condition)` that returns a named numeric vector of
#'   parameter values.
#' @param simulate `function(truth, condition)` that returns the data.
#'
#' @return A generator. Each dataset is `list(truth = , data = )`.
#'   Built-in measurements read `truth` from it.
#' @export
prior_draws <- function(prior, simulate) {
  .sg_function(prior, "prior")
  .sg_function(simulate, "simulate")
  function(condition, context) {
    truth <- prior(condition)
    if (!is.numeric(truth) || is.null(names(truth)) || anyNA(names(truth))) {
      stop("`prior` must return a named numeric vector.")
    }
    list(truth = truth, data = simulate(truth, condition))
  }
}

#' A conjugate normal linear model
#'
#' `conjugate_lm()` is a [study_method()] for the model
#' `y ~ normal(X beta, sigma)` with a known `sigma` and independent
#' `normal(prior_mean, prior_sd)` priors on `beta`. The posterior is exact,
#' so the method needs no sampler. It is useful to test a study design and to
#' check calibration.
#'
#' The method expects datasets of the form that [prior_draws()] makes, where
#' `data` is `list(X = , y = )`. The column names of `X` are the parameter
#' names.
#'
#' @param prior_mean,prior_sd Prior mean and standard deviation of each
#'   coefficient.
#' @param sigma Known residual standard deviation.
#' @param draws Number of posterior draws for the `draws` artifact.
#'
#' @return A `bayesim_study_method` with the artifact `draws`, a matrix with
#'   one column per coefficient.
#' @export
conjugate_lm <- function(prior_mean = 0, prior_sd = 1, sigma = 1, draws = 1000) {
  study_method(
    fit = function(data, context) {
      s <- context$settings
      X <- data$data$X
      y <- data$data$y
      if (is.null(colnames(X))) {
        stop("`X` needs column names, the parameter names.")
      }
      precision <- crossprod(X) / s$sigma^2 + diag(1 / s$prior_sd^2, ncol(X))
      cov <- chol2inv(chol(precision))
      mean <- drop(cov %*% (crossprod(X, y) / s$sigma^2 + s$prior_mean / s$prior_sd^2))
      dimnames(cov) <- list(colnames(X), colnames(X))
      list(mean = stats::setNames(mean, colnames(X)), cov = cov)
    },
    extract = list(
      draws = function(fit, data, context) {
        p <- length(fit$mean)
        n <- context$settings$draws
        z <- matrix(stats::rnorm(n * p), n, p)
        out <- z %*% chol(fit$cov) + rep(fit$mean, each = n)
        colnames(out) <- names(fit$mean)
        out
      }
    ),
    settings = list(
      prior_mean = prior_mean,
      prior_sd = prior_sd,
      sigma = sigma,
      draws = draws
    )
  )
}

# Built-in measurements ----------------------------------------------------------------

# Draws as a matrix with one column per variable. Accepts a matrix, an array
# of iterations x chains x variables, or a posterior draws object.
.draws_matrix <- function(draws) {
  if (inherits(draws, "draws")) {
    draws <- posterior::as_draws_matrix(draws)
    keep <- !colnames(draws) %in% c(".chain", ".iteration", ".draw")
    return(unclass(draws)[, keep, drop = FALSE])
  }
  if (is.array(draws) && length(dim(draws)) == 3L) {
    d <- dim(draws)
    out <- matrix(draws, d[1] * d[2], d[3])
    colnames(out) <- dimnames(draws)[[3]]
    return(out)
  }
  if (is.matrix(draws)) {
    return(draws)
  }
  stop("Draws must be a matrix, a 3-D array or a posterior draws object.")
}

# Draws and truth for the parameters that both name.
.draws_truth <- function(artifacts) {
  draws <- .draws_matrix(artifacts$draws)
  truth <- artifacts$data$truth
  if (is.null(truth) || is.null(names(truth))) {
    stop("The dataset needs a named `truth`, as `prior_draws()` makes it.")
  }
  if (is.null(colnames(draws))) {
    stop("The draws need variable names.")
  }
  target <- intersect(colnames(draws), names(truth))
  if (!length(target)) {
    stop("The draws and the truth share no parameter names.")
  }
  list(draws = draws[, target, drop = FALSE], truth = truth[target], target = target)
}

.builtin_measure <- function(f, summary = "mean") {
  structure(f, needs = c("data", "draws"), summary = summary)
}

#' Built-in measurements
#'
#' Measurements of parameter recovery for [with_measure()]. Each one needs
#' the artifacts `draws` (posterior draws, with one column per parameter) and
#' `data` (a dataset with a named `truth`, as [prior_draws()] makes it). It
#' returns one row per parameter, with columns `target` and `value`.
#'
#' * `measure_coverage(level)`: `value` is 1 when the central posterior interval
#'   holds the truth, otherwise 0.
#' * `measure_bias()`: `value` is the posterior mean minus the truth.
#' * `measure_rmse()`: `value` is the squared error of the posterior mean.
#'   [summarise_study()] reports the root of its mean.
#' * `measure_sbc_rank()`: `value` is the rank of the truth among the draws divided
#'   by the number of draws. When the method is calibrated, it is uniform on
#'   0 to 1 and its mean is 0.5. The columns `rank` and `draws` hold the rank
#'   and the number of draws.
#'
#' @param level Interval probability.
#' @return A measurement function with default `needs`.
#' @name study-measures
NULL

#' @rdname study-measures
#' @export
measure_coverage <- function(level = 0.9) {
  if (!is.numeric(level) || length(level) != 1L || !(level > 0 && level < 1)) {
    cli::cli_abort("{.arg level} must be a probability between 0 and 1.")
  }
  .builtin_measure(function(artifacts, context) {
    x <- .draws_truth(artifacts)
    p <- (1 - level) / 2
    lower <- apply(x$draws, 2L, stats::quantile, probs = p, names = FALSE)
    upper <- apply(x$draws, 2L, stats::quantile, probs = 1 - p, names = FALSE)
    data.frame(
      target = x$target,
      value = as.numeric(lower <= x$truth & x$truth <= upper),
      row.names = NULL
    )
  })
}

#' @rdname study-measures
#' @export
measure_bias <- function() {
  .builtin_measure(function(artifacts, context) {
    x <- .draws_truth(artifacts)
    data.frame(
      target = x$target,
      value = colMeans(x$draws) - x$truth,
      row.names = NULL
    )
  })
}

#' @rdname study-measures
#' @export
measure_rmse <- function() {
  .builtin_measure(
    function(artifacts, context) {
      x <- .draws_truth(artifacts)
      data.frame(
        target = x$target,
        value = (colMeans(x$draws) - x$truth)^2,
        row.names = NULL
      )
    },
    summary = "rmse"
  )
}

#' @rdname study-measures
#' @export
measure_sbc_rank <- function() {
  .builtin_measure(function(artifacts, context) {
    x <- .draws_truth(artifacts)
    n <- nrow(x$draws)
    rank <- colSums(sweep(x$draws, 2L, x$truth, `<`))
    data.frame(
      target = x$target,
      value = rank / n,
      rank = as.integer(rank),
      draws = n,
      row.names = NULL
    )
  })
}

#' Paired difference between methods
#'
#' `compare_paired()` makes a comparison for [with_comparison()]. Within
#' each group, usually one dataset, it subtracts the reference method's
#' measurement from each other method's, target by target. Because both
#' methods saw the same data, the differences have a smaller Monte Carlo
#' error than the difference of two separate means. [summarise_study()]
#' reports that paired error.
#'
#' @param measurement Name of the measurement to compare.
#' @param reference Method ID of the reference method.
#'
#' @return A comparison function. Its rows have the columns `method_id`,
#'   `reference`, `target` and `value`. A dataset where the reference failed
#'   gives no rows.
#' @export
compare_paired <- function(measurement, reference) {
  .sg_name(measurement, "measurement")
  .sg_name(reference, "reference")
  function(results, artifacts, context) {
    rows <- results[results$measurement == measurement, , drop = FALSE]
    if (!"target" %in% names(rows)) {
      rows$target <- rep(NA_character_, nrow(rows))
    }
    ref <- rows[rows$method_id == reference, , drop = FALSE]
    other <- rows[rows$method_id != reference, , drop = FALSE]
    hit <- match(other$target, ref$target)
    other <- other[!is.na(hit), , drop = FALSE]
    hit <- hit[!is.na(hit)]
    data.frame(
      method_id = other$method_id,
      reference = rep(reference, nrow(other)),
      target = other$target,
      value = other$value - ref$value[hit],
      stringsAsFactors = FALSE
    )
  }
}
