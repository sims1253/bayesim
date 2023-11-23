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

#' Full forward sampling of a the response of brms fit, including multivariate models.
#'
#' @source This function is taken from the SBC package
#'   (https://github.com/hyunjimoon/SBC) and only here to ensure stability, as
#'   it is not exported in SBC.
#'
#' @param fit An object of class `brmsfit`
#' @param newdata An optional data.frame for which to evaluate predictions. If NULL (default), the original data of the model is used.
#' @param draws An integer vector specifying the posterior draws to be used. If NULL (the default), all draws are used.
#' @param validate_all if TRUE, validation of input data will be done in all iterations, otherwise only once
#'
#' @return A list of data.frames containing the draws.
#'
#' @keywords internal
#' @examples # Pending
#' @export
brms_full_ppred <- function(fit, newdata = NULL, draws = NULL, validate_all = FALSE) {
  # 1. determine term hierarchy
  resp <- brms_response_sequence(fit)
  # 2.1. initialize dataframe using the original fit's data
  if (is.null(newdata)) newdata <- fit$data
  n <- nrow(newdata)
  # 2.3. if no draws set, range from 1 to all iters (check draws < iters)
  if (is.null(draws)) draws <- seq_len(sum(fit$fit@sim$n_save))
  # 2.4. create list to hold data
  pp_data <- list()

  if (!validate_all) {
    # Validate once
    newdata <- brms::validate_newdata(newdata, fit)
  }

  for (i in draws) {
    pp_data[[i]] <- newdata
    for (vars in resp) {
      pp_data[[i]][, vars] <- array(
        brms::posterior_predict(
          fit,
          newdata = pp_data[[i]],
          resp = vars, draw_ids = i,
          skip_validate = !validate_all
        ),
        dim = c(1, n, length(vars))
      )[1, , ]
    }
  }
  pp_data
}

nodes_by_depth <- function(adj_matrix) {
  depth_list <- list()
  var_names <- rownames(adj_matrix)
  while (nrow(adj_matrix)) {
    pos <- which(apply(adj_matrix, 1, sum) == 0)
    depth_list <- c(depth_list, list(var_names[pos]))

    var_names <- var_names[-pos]
    adj_matrix <- adj_matrix[-pos, -pos, drop = FALSE]
  }
  depth_list
}

#' Determine the response sequence of brms model
#' @source This function is taken from the SBC package
#'   (https://github.com/hyunjimoon/SBC) and only here to ensure stability, as
#'   it is not exported in SBC.
#' @keywords internal
brms_response_sequence <- function(x) {
  UseMethod("brms_response_sequence")
}

#' @source This function is taken from the SBC package
#'   (https://github.com/hyunjimoon/SBC) and only here to ensure stability, as
#'   it is not exported in SBC.
#' @method brms_response_sequence brmsfit
#' @export
brms_response_sequence.brmsfit <- function(x) {
  brms_response_sequence(x$formula)
}

#' @source This function is taken from the SBC package
#'   (https://github.com/hyunjimoon/SBC) and only here to ensure stability, as
#'   it is not exported in SBC.
#' @method brms_response_sequence bform
#' @export
brms_response_sequence.bform <- function(x) {
  term_list <- brms_response_sequence(brms::brmsterms(x))
  resp_vars <- names(term_list)

  adjacency <- t(sapply(term_list, \(x)is.element(resp_vars, x)))
  attr(adjacency, "dimnames") <- list(resp_vars, resp_vars)
  nodes_by_depth(adjacency)
}

#' @source This function is taken from the SBC package
#'   (https://github.com/hyunjimoon/SBC) and only here to ensure stability, as
#'   it is not exported in SBC.
#' @method brms_response_sequence mvbrmsterms
#' @export
brms_response_sequence.mvbrmsterms <- function(x) {
  names(x$terms) <- NULL
  sapply(x$terms, brms_response_sequence)
}

#' @source This function is taken from the SBC package
#'   (https://github.com/hyunjimoon/SBC) and only here to ensure stability, as
#'   it is not exported in SBC.
#' @method brms_response_sequence brmsterms
#' @export
brms_response_sequence.brmsterms <- function(x) {
  vars <- list(unique(unlist(lapply(x$dpars, brms_response_sequence))))
  names(vars) <- all.vars(x$respform)
  vars
}

#' @source This function is taken from the SBC package
#'   (https://github.com/hyunjimoon/SBC) and only here to ensure stability, as
#'   it is not exported in SBC.
#' @method brms_response_sequence btl
#' @export
brms_response_sequence.btl <- function(x) {
  c("1", all.vars(x$formula))
}
