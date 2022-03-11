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
      "elpd" = loo_object$estimates[1, 1],
      "se_elpd" = loo_object$estimates[1, 2],
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
  final_result <- data.frame(matrix(nrow = length(loo_object_matrix), ncol = 2*length(predictive_metrics)))
  colnames(final_result) <- unlist(lapply(predictive_metrics, function(x){c(paste0(x, "_delta"), paste0(x, "_se_delta"))}))
  index <- names(loo_object_matrix)
  rownames(final_result) <- index

  for (i in seq_along(predictive_metrics)) {
    metric <- predictive_metrics[[i]]

    loo_result <- loo_compare(lapply(loo_object_matrix, function(x){x[[i]]}))
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

