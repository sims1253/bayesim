#' Simulate a new dataset using forward sampling.
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
forward_sampling <- function(fit, i, n) {
  if (i > ceiling(
    (fit$fit@sim$iter - fit$fit@sim$warmup) / fit$fit@sim$thin * fit$fit@sim$chains
  ) | i <= 0) {
    stop("You tried to use a non-existent posterior sample.")
  }

  model_formula <- fit$formula
  if (is(fit$formula, "brmsformula")) {
    response <- model_formula$resp
    if (length(all.vars(model_formula$formula)) > 1) {
      stop("Your model does not contain responses for all used variables.")
    }
    df <- data.frame(matrix(0, # Create empty data frame
      nrow = n,
      ncol = 1
    ))
    colnames(df) <- c(response)
    df[response] <- as.vector(
      brms::posterior_predict(fit,
        newdata = df[response],
        resp = response,
        draw_ids = i
      )
    )
  } else {
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
  }
  return(df)
}
