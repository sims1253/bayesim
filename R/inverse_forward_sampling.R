#' Simulate a new dataset using forward sampling.
#'
#' @param fit An object of class brmsfit or a list of brmsfit objects.
#' @param i The index of a single posterior draw to simulate a dataset for.
#'  The index is passed to \code{\link{posterior_predict}}'s "draw_ids"
#'  argument.
#' @param newdata A dataframe that is passed to posterior predict.
#' @param ... Potential additional arguments.
#'
#' @return A data.frame containing n observations for each variable in the fit.
#' @export forward_sampling forward_sampling.list
#'
#'
forward_sampling <- function(fit, i, n, ...) {
  UseMethod("forward_sampling")
}

#'
forward_sampling.list <- function(x, i, n, ...) {
  if (any(lapply(x, bayeshear::post_warmup_samples) < i) | i <= 0) {
    stop("You tried to use a non-existent posterior sample.")
  }

  resp_list <- list()
  var_list <- list()

  for (index in seq_along(x)) {
    fit <- x[[index]]
    if (is(fit$formula, "brmsformula")) {
      resp_list[[fit$formula$resp]] <- index
      var_list[[fit$formula$resp]] <- all.vars(fit$formula$formula)[-1]
    } else if (is(fit$formula, "mvbrmsformula")) {
      for (formula in fit$formula$forms) {
        resp_list[[formula$resp]] <- index
        var_list[[formula$resp]] <- all.vars(formula$formula)[-1]
      }
    } else {
      stop("Unsupported model type detected!
           Please use brmsfit or stanfit objects.")
    }
  }
  df <- data.frame(matrix(0, # Create empty data frame
    nrow = n,
    ncol = length(resp_list)
  ))
  colnames(df) <- names(resp_list)

  missing_responses <- names(resp_list)
  while (length(missing_responses) > 0) {
    success <- FALSE
    for (response in missing_responses) {
      if (length(intersect(missing_responses, var_list[[response]])) == 0) {
        df[response] <- as.vector(
          brms::posterior_predict(x[[resp_list[[response]]]],
            newdata = df[var_list[[response]]],
            resp = response,
            draw_ids = i,
            ...
          )
        )
        missing_responses <- missing_responses[!(missing_responses == response)]
        success <- TRUE
      }
    }
    if (!success) {
      stop("Unsolvable variable dependencies detected.")
    }
  }
  return(df)
}
