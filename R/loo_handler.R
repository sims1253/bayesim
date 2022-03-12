#' Title
#'
#' @param fit
#'
#' @return
#' @export
#'
#' @examples
elpd_loo_handler <- function(fit) {
  loo_object <- brms::loo(fit, save_psis = TRUE)
  return(
    list(
      "p_loo" = loo_object$estimates[2, 1],
      "se_p_loo" = loo_object$estimates[2, 2],
      "elpd_loo" = loo_object$estimates[1, 1],
      "se_elpd_loo" = loo_object$estimates[1, 2],
      "looic" = loo_object$estimates[3, 1],
      "se_looic" = loo_object$estimates[3, 2],
      "object" = loo_object
    )
  )
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
  final_result <- data.frame(matrix(nrow = length(loo_object_matrix), ncol = 2 * length(predictive_metrics)))
  colnames(final_result) <- unlist(lapply(predictive_metrics, function(x) {
    c(paste0(x, "_delta"), paste0(x, "_se_delta"))
  }))
  index <- names(loo_object_matrix)
  rownames(final_result) <- index

  for (i in seq_along(predictive_metrics)) {
    metric <- predictive_metrics[[i]]

    loo_result <- loo_compare(lapply(loo_object_matrix, function(x) {
      x[[i]]
    }))
    deltas <- numeric(length = length(index))
    errors <- numeric(length = length(index))
    for (i in seq_along(index)) {
      deltas[[i]] <- loo_result[index[[i]], "elpd_diff"]
      errors[[i]] <- loo_result[index[[i]], "se_diff"]
    }
    final_result[paste0(metric, "_delta")] <- deltas
    final_result[paste0(metric, "_se_delta")] <- errors
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
#' @param psis_object \code{brms:::.psis} object for psis diagnostics
#'
#' @return a loo object, containing a criterion, disguised as elpd
#' @export
#'
#' @examples
custom_loo_object <- function(pointwise_criterion,
                              psis_object = NULL) {
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
#' @return \code{custom_loo_object} object with rmse acting as elpd.
#' @export
#'
#' @examples
rmse_loo <- function(fit,
                     psis_object = NULL,
                     ...) {
  if (is.null(psis_object)) {
    psis_object <- brms:::.psis(fit, newdata = fit$data, resp = NULL)
  }
  pointwise_rmse <- rmse(
    y = brms::get_y(fit),
    yrep = posterior_predict(fit, fit$data),
    weights = exp(psis_object$log_weights)
  )

  loo_object <- custom_loo_object(
    pointwise_criterion = -pointwise_rmse,
    psis_object = psis_object
  )
  return(
    list(
      "rmse_loo" = loo_object$estimates[1, 1],
      "se_rmse_loo" = loo_object$estimates[1, 2],
      "object" = loo_object
    )
  )
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
rmse_newdata <- function(fit, newdata) {
  pointwise_rmse <- rmse(
    y = y <- newdata$y,
    yrep = brms::posterior_predict(fit, newdata = newdata)
  )

  loo_object <- custom_loo_object(
    pointwise_criterion = -pointwise_rmse
  )
  return(
    list(
      "rmse_newdata" = loo_object$estimates[1, 1],
      "se_rmse_newdata" = loo_object$estimates[1, 2],
      "object" = loo_object
    )
  )
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
    return(sqrt(colMeans(((y_matrix - yrep)^2
    ))))
  } else {
    return(sqrt(colSums(weights * ((y_matrix - yrep)^2
    )) /
      colSums(weights)))
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
elpd_newdata <- function(fit, newdata) {
  ll <- brms::log_lik(fit, newdata = newdata)
  elpd <- matrixStats::colLogSumExps(ll) - log(nrow(ll))
  loo_object <- custom_loo_object(
    pointwise_criterion = elpd
  )
  return(
    list(
      "elpd_newdata" = loo_object$estimates[1, 1],
      "se_elpd_newdata" = loo_object$estimates[1, 2],
      "object" = loo_object
    )
  )
}
