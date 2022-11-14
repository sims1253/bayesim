#' Title
#'
#' @param numeric_metrics
#' @param predictive_metrics
#' @param testing_data...
#' @param fit
#'
#' @return
#' @export
#'
#' @examples
metric_list_handler <- function(fit,
                                numeric_metrics,
                                predictive_metrics,
                                testing_data,
                                ...) {
  posterior_draws <- posterior::extract_variable_matrix(fit, variable = "b_x")

  loo_objects <- vector(mode = "list", length = length(predictive_metrics))
  predictive_results <- vector(mode = "list", length = length(predictive_metrics))

  psis_object <- NULL
  for (i in seq_along(predictive_metrics)) {
    identifier <- predictive_metrics[[i]]
    result <- metric_lookup(
      identifier = identifier,
      fit = fit,
      posterior_draws = posterior_draws,
      testing_data = testing_data,
      psis_object,
      ...
    )
    if (is.null(psis_object) & "psis_object" %in% names(result$object)) {
      psis_object <- result$object$psis_object
    }
    loo_objects[[i]] <- result$object
    predictive_results[[i]] <- result[names(result) != "object"]
  }

  numeric_result_list <- vector(mode = "list", length = length(numeric_metrics))
  for (i in seq_along(numeric_metrics)) {
    identifier <- numeric_metrics[[i]]
    result <- metric_lookup(
      identifier = identifier,
      fit = fit,
      posterior_draws = posterior_draws,
      psis_object = psis_object,
      ...
    )
    if (length(result) == 1) {
      result <- list(result)
      names(result) <- identifier
      numeric_result_list[[i]] <- result
    } else {
      numeric_result_list[[i]] <- result
    }
  }

  return(
    list(
      numeric_results = unlist(c(numeric_result_list, predictive_results),
        recursive = FALSE
      ),
      loo_objects = loo_objects
    )
  )
}
