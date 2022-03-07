#' Title
#'
#' @param fit
#' @param metric_list
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
metric_list_handler <- function(fit, metric_list, ...) {
  posterior_draws <- as.vector(posterior::as_draws_array(fit, variable = "b_x"))
  result_list <- list()
  for (i in seq_along(metric_list)) {
    identifier <- metric_list[[i]]
    result <- metric_lookup(
      identifier = identifier,
      fit = fit,
      posterior_draws = posterior_draws,
      ...
    )

    if (length(result) == 1) {
      list_build <- list(result)
      names(list_build) <- paste(identifier)
      result_list <- append(result_list, list_build)
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
        result_list <- append(result_list, list_build)
      }
    }
  }
  return(result_list)
}
