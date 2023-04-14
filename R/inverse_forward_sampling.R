#' Simulate a new dataset using inverse forward sampling.
#'
#' @param fit An object of class brmsfit.
#' @param i The index of a single posterior draw to simulate a dataset for.
#'          The index is passed to posterior_predict's "draw_ids" argument.
#' @param n The number of samples for the newly simulated dataset.
#'
#' @return A data.frame containing n observations for each variable in the fit.
#' @export
#'
#' @examples
inverse_forward_sampling <- function(fit, i, n) {
  model_formula <- fit$formula
  df <- data.frame(matrix(0, # Create empty data frame
    nrow = n,
    ncol = length(model_formula$responses)
  ))
  colnames(df) <- unname(model_formula$responses)
  missing_responses <- unname(model_formula$responses)
  while (length(missing_responses) > 0) {
    for (response in missing_responses) {
      cur_formula <- model_formula$forms[response]
      all_variables <- all.vars(cur_formula[[response]]$formula)
      if (length(intersect(missing_responses, all_variables)) == 1) {
        missing_variables <- all_variables[
          !(all_variables %in% unname(model_formula$responses))
        ]
        if (length(missing_variables) != 0) {
          stop(
            paste(
              missing_variables,
              "are needed to model",
              response,
              "but are not modeled by the provided model."
            )
          )
        }
        newdata <- df[all_variables[!(all_variables == response)]]
        df[response] <- as.vector(
          brms::posterior_predict(fit,
            newdata = newdata,
            resp = response,
            draw_ids = i
          )
        )
        missing_responses <- missing_responses[!(missing_responses == response)]
      }
    }
  }
  return(df)
}
