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
  predictive_results <- vector(mode = "list", length = length(numeric_metrics))
  numeric_results <- vector(mode = "list", length = length(numeric_metrics))

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
    if ("psis_object" %in% names(result$object)){
      psis_object <- result$object$psis_object
    }
    loo_objects[[i]] <-result$object
    predictive_results[[i]] <- result[names(result) != "object"]
  }

  numeric_results <- vector(mode = "list", length = length(numeric_metrics))
  for (i in seq_along(numeric_metrics)) {
    identifier <- numeric_metrics[[i]]
    result <- metric_lookup(
      identifier = identifier,
      fit = fit,
      posterior_draws = posterior_draws,
      predictive_results = predictive_results,
      ...
    )

    if (length(result) == 1) {
      list_build <- list(result)
      names(list_build) <- paste(identifier)
      numeric_results[[i]] <- list_build
    } else {
      for (j in seq_along(result)) {
        list_build <- list(result[[names(result)[[j]]]])
        if (grepl("<|>|=", identifier)) {
          identifier <- gsub(">", "_gr_", identifier)
          identifier <- gsub("<", "_sm_", identifier)
          identifier <- gsub("=", "_eq_", identifier)
        }
        names(list_build) <- paste(identifier,
          names(result)[[j]],
          sep = "_"
        )
        numeric_results[[i]] <- list_build
      }
    }
  }


  return(
    list(
      numeric_results = unlist(c(numeric_results, predictive_results)),
      loo_objects = loo_objects
    )
  )
}
