#' Title
#'
#' @param fit
#'
#' @return
#' @export
#'
#' @examples
elpd_loo_handler <- function(fit) {
  loo_object <- brms::loo(fit)
  return(
    list(
      "p_loo" = loo_object$estimates[2, 1],
      "se_p_loo" = loo_object$estimates[2, 2],
      "elpd_loo" = loo_object$estimates[1, 1],
      "se_elpd_loo" = loo_object$estimates[1, 2],
      "looic" = loo_object$estimates[3, 1],
      "se_looic" = loo_object$estimates[3, 2]
    )
  )
}

#' Pointwise elps summaries
#'
#' Convenience function to collect quantiles and summaries of the pointwise elpd
#' estimates instead of just the main estimates.
#' The returned summaries are all sample size independent.
#'
#' @param fit A brmsfit object.
#' @param quantiles A vector of quantiles of interest.
#' @param newdata If supplied, returns the summaries for [elpd_test()]
#'                otherwise, returns [brms::elpd()] summaries by
#'                default.
#' @return A named list of summaries.
#' @export
#'
#' @examples
#' fit <- brms::brm(y ~ 1, data = list(rnorm(1000)))
#' elpd_pointwise_summaries(fit, seq(0.1, 0.9, length.out = 9))
#' elpd_pointwise_summaries(fit, seq(0.1, 0.9, length.out = 9), rnorm(1000))
elpd_pointwise_summaries <- function(fit, quantiles, newdata = NULL) {
  if (is.null(newdata)) {
    loo_object <- brms::loo(fit)

    elpd_pointwise <- loo_object$pointwise[, 1]
    elpd_quantiles <- as.list(quantile(elpd_pointwise, prob = quantiles))
    names(elpd_quantiles) <- lapply(
      quantiles,
      function(x) paste0("elpd_loo_quantile_", x)
    )

    p_loo_pointwise <- loo_object$pointwise[, 3]
    p_loo_quantiles <- as.list(quantile(p_loo_pointwise, prob = quantiles))
    names(p_loo_quantiles) <- lapply(
      quantiles,
      function(x) paste0("p_loo_quantile_", x)
    )
    out <- c(elpd_quantiles, p_loo_quantiles)
    out$elpd_loo_mean <- mean(elpd_pointwise)
    out$elpd_loo_sd <- sd(elpd_pointwise)
    out$p_loo_mean <- sd(p_loo_pointwise)
    out$p_loo_sd <- sd(p_loo_pointwise)
    out$elpd_loo_se_mean <- mean(loo_object$pointwise[, 2])
    return(out)
  } else {
    loo_object <- elpd_test(fit, newdata, TRUE)

    elpd_pointwise <- loo_object$pointwise[, 1]
    elpd_quantiles <- as.list(quantile(elpd_pointwise, prob = quantiles))
    names(elpd_quantiles) <- lapply(
      quantiles,
      function(x) paste0("elpd_test_quantile_", x)
    )
    out <- c(elpd_quantiles)
    out$elpd_test_mean <- mean(elpd_pointwise)
    out$elpd_test_sd <- sd(elpd_pointwise)
    out$elpd_test_se_mean <- loo_object$estimates[1, 2] / length(elpd_pointwise)
    return(out)
  }
}

#' Title
#'
#' @param loo_object_matrix
#' @param predictive_metrics
#'
#' @return
#' @export
#'
#' @examples
loo_compare_handler <- function(loo_object_matrix, predictive_metrics) {
  compare_metrics <- substr(
    predictive_metrics[
      grep(
        "_compare",
        predictive_metrics
      )
    ],
    1,
    nchar(predictive_metrics[1]) - 8
  )

  final_result <- data.frame(matrix(
    nrow = length(loo_object_matrix),
    ncol = 2 * length(predictive_metrics)
  ))
  colnames(final_result) <- unlist(lapply(compare_metrics, function(x) {
    c(paste0(x, "_delta"), paste0(x, "_se_delta"))
  }))
  index <- names(loo_object_matrix)
  rownames(final_result) <- index

  valid_entries <- c()
  for (i in seq_along(loo_object_matrix)) {
    if (
      !is.null(loo_object_matrix[[i]]) &
        !any(sapply(loo_object_matrix[[i]], is.null))
    ) {
      valid_entries <- c(valid_entries, i)
    }
  }

  if (length(valid_entries) < 2) {
    for (i in seq_along(compare_metrics)) {
      metric <- compare_metrics[[i]]
      deltas <- numeric(length = length(index))
      errors <- numeric(length = length(index))
      for (i in seq_along(index)) {
        deltas[[i]] <- NA
        errors[[i]] <- NA
      }
      final_result[paste0(metric, "_delta")] <- deltas
      final_result[paste0(metric, "_se_delta")] <- errors
    }
  } else {
    for (i in seq_along(compare_metrics)) {
      metric <- compare_metrics[[i]]

      loo_result <- brms::loo_compare(lapply(
        loo_object_matrix[valid_entries],
        function(x) {
          x[[i]]
        }
      ))
      deltas <- numeric(length = length(index))
      errors <- numeric(length = length(index))

      for (i in seq_along(index)) {
        if (
          any(sapply(loo_object_matrix[[i]], is.null)) |
            is.null(loo_object_matrix[[i]])
        ) {
          deltas[[i]] <- NA
          errors[[i]] <- NA
        } else {
          deltas[[i]] <- loo_result[index[[i]], "elpd_diff"]
          errors[[i]] <- loo_result[index[[i]], "se_diff"]
        }
      }
      final_result[paste0(metric, "_delta")] <- deltas
      final_result[paste0(metric, "_se_delta")] <- errors
    }
  }
  return(final_result)
}


#' Builds a loo object that contains any pointwise criterion, acting as elpd
#' for compatibility.
#'
#' LOO currently has hardcoded expectations of elpd as part of loo objects
#' so to use loo objects, we have to disguise other criterions as elpd.
#'
#' @param pointwise_criterion vector of criterion values for each observation
#' @param psis_object `brms:::.psis` object for psis diagnostics
#'
#' @return a loo object, containing a criterion, disguised as elpd
#' @export
#'
#' @examples
custom_loo_object <- function(pointwise_criterion, psis_object = NULL) {
  loo_object <- list()
  criterion <- sum(pointwise_criterion)
  se_criterion <- sqrt(length(pointwise_criterion) * var(pointwise_criterion))
  loo_object$estimates <- matrix(c(criterion, se_criterion), nrow = 1)
  loo_object$pointwise <- as.matrix(pointwise_criterion)
  loo_object$elpd_loo <- criterion

  row.names(loo_object$estimates) <- c("elpd_loo")
  colnames(loo_object$estimates) <- c("Estimate", "SE")
  colnames(loo_object$pointwise) <- c("elpd_loo")

  if (!is.null(psis_object)) {
    loo_object$diagnostics <- psis_object$diagnostics
  }

  attr(loo_object, "model_name") <- NULL
  attr(loo_object, "dims") <- dim(psis_object)
  attr(loo_object, "class") <-
    c("psis_loo", "importance_sampling_loo", "loo")

  return(loo_object)
}

#' Calculate PSIS-loo rmse for a given brms fit
#'
#' @param ... Additional arguments to be passed to update() in case of reloo
#' @param fit
#' @param psis_object
#'
#' @return `custom_loo_object` object with rmse acting as elpd.
#' @export
#'
#' @examples
rmse_loo <- function(
  fit,
  psis_object = NULL,
  return_object = FALSE,
  yrep = NULL,
  ...
) {
  if (is.null(psis_object)) {
    psis_object <- brms:::.psis(fit, newdata = fit$data, resp = NULL)
  }
  if (is.null(yrep)) {
    yrep <- brms::posterior_predict(fit, fit$data)
  }
  pointwise_rmse <- rmse(
    y = brms::get_y(fit),
    yrep = yrep,
    weights = exp(psis_object$log_weights)
  )
  if (return_object) {
    return(custom_loo_object(pointwise_criterion = pointwise_rmse))
  } else {
    return(
      list(
        "rmse_loo" = sum(pointwise_rmse),
        "se_rmse_loo" = sqrt(length(pointwise_rmse) * var(pointwise_rmse))
      )
    )
  }
}

#' Title
#'
#' @param y
#' @param yrep
#'
#' @return
#' @export
#'
#' @examples
rmse_test <- function(fit, newdata, return_object = FALSE) {
  pointwise_rmse <- rmse(
    y = y <- newdata$y,
    yrep = brms::posterior_predict(fit, newdata = newdata)
  )

  if (return_object) {
    return(custom_loo_object(pointwise_criterion = pointwise_rmse))
  } else {
    return(
      list(
        "rmse_test" = sum(pointwise_rmse),
        "se_rmse_test" = sqrt(length(pointwise_rmse) * var(pointwise_rmse))
      )
    )
  }
}

#' Calculate the root-mean-squared-error for given y and yrep.
#'
#' If psis-weights are supplied, the corresponding psis-rmse is returned.
#'
#' @param y Vector of observed values
#' @param yrep Vector of predicted Values
#' @param weights PSIS weights
#'
#' @return rmse for the given y and yrep vectors
#' @export
#'
#' @examples
rmse <- function(y, yrep, weights = NULL) {
  y_matrix <- matrix(
    y,
    nrow <- nrow(yrep),
    ncol <- ncol(yrep),
    byrow <- TRUE
  )
  if (is.null(weights)) {
    return(sqrt(colMeans(((y_matrix - yrep)^2))))
  } else {
    return(sqrt(
      colSums(weights * ((y_matrix - yrep)^2)) /
        colSums(weights)
    ))
  }
}

#' Title
#'
#' @param fit
#' @param newdata
#'
#' @return
#' @export
#'
#' @examples
elpd_test <- function(fit, newdata, return_object = FALSE) {
  ll <- brms::log_lik(fit, newdata = newdata)
  elpd <- matrixStats::colLogSumExps(ll) - log(nrow(ll))
  if (return_object) {
    return(custom_loo_object(pointwise_criterion = elpd))
  } else {
    return(
      list(
        "elpd_test" = sum(elpd),
        "se_elpd_test" = sqrt(length(elpd) * var(elpd))
      )
    )
  }
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
r2 <- function(y, yrep, weights = NULL) {
  ss_y <- sum((y - mean(y))^2)
  pointwise_loo_r2 <- vector(mode = "numeric", length = length(y))

  if (is.null(weights)) {
    for (n in seq_along(pointwise_loo_r2)) {
      ss_e <- sum((y[n] - yrep[, n])^2)
      pointwise_loo_r2[[n]] <- 1 / length(y) - ss_e / ss_y
    }
  } else {
    for (n in seq_along(pointwise_loo_r2)) {
      ss_e <- sum((weights[, n] * (y[n] - yrep[, n])^2) / sum(weights[, n]))
      pointwise_loo_r2[[n]] <- 1 / length(y) - ss_e / ss_y
    }
  }
  return(pointwise_loo_r2)
}

#' Calculate PSIS-loo R² for a given brms fit
#'
#' @param fit brms fit to calculate rmse for
#' @param psis_object
#' @param ...
#'
#' @return `custom_loo_object` object with R² acting as elpd.
#' @export
#'
#' @examples
r2_loo <- function(
  fit,
  psis_object = NULL,
  yrep = NULL,
  return_object = FALSE,
  ...
) {
  if (is.null(psis_object)) {
    psis_object <- brms:::.psis(fit, newdata = fit$data, resp = NULL)
  }
  if (is.null(yrep)) {
    yrep <- brms::posterior_predict(fit, fit$data)
  }
  pointwise_loo_r2 <- r2(
    y = brms::get_y(fit),
    yrep = yrep,
    weights = exp(psis_object$log_weights)
  )
  if (return_object) {
    return(custom_loo_object(pointwise_criterion = pointwise_loo_r2))
  } else {
    return(
      list(
        "r2_loo" = sum(pointwise_loo_r2),
        "se_r2_loo" = sqrt(length(pointwise_loo_r2) * var(pointwise_loo_r2))
      )
    )
  }
}

#' Title
#'
#' @param fit
#' @param newdata
#'
#' @return
#' @export
#'
#' @examples
r2_test <- function(fit, newdata, return_object = FALSE) {
  pointwise_r2 <- r2(
    y = y <- newdata$y,
    yrep = brms::posterior_predict(fit, newdata = newdata)
  )

  if (return_object) {
    return(custom_loo_object(pointwise_criterion = pointwise_r2))
  } else {
    return(
      list(
        "r2_test" = sum(pointwise_r2),
        "se_r2_test" = sqrt(length(pointwise_r2) * var(pointwise_r2))
      )
    )
  }
}
