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
  result_list <- list()
  for (i in seq_along(metric_list)) {
    identifier <- metric_list[[i]]
    result <- metric_lookup(
      identifier = identifier,
      fit = fit,
      ...
    )

    if (length(result) == 1) {
      list_build <- list(result)
      names(list_build) <- paste(identifier)
      result_list <- append(result_list, list_build)
    } else {
      for (j in seq_along(result)) {
        list_build <- list(result[[names(result)[[j]]]])
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
