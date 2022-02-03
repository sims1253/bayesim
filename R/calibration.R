#' Title
#'
#' @param model_name
#' @param results
#' @param samples
#'
#' @return
#' @export
#'
#' @examples
calibration_result_text <- function(model_name,
                                    results,
                                    samples,
                                    hypothesis_list) {
  print(model_name)

  for (i in seq_along(hypothesis_list)) {
    current <- hypothesis_list[[i]]

    lower_ci <- unlist(
      lapply(
        results[[current]],
        function(x) {
          x$hypothesis$CI.Lower
        }
      )
    )
    upper_ci <- unlist(
      lapply(
        results[[current]],
        function(x) {
          x$hypothesis$CI.Upper
        }
      )
    )

    if (current == "x=0") {
      print(
        paste(
          sum(lower_ci > 0 | upper_ci < 0),
          "/",
          samples,
          " of 'x=0' 95% CIs exclude 0.",
          sep = ""
        )
      )
    }

    if (current == "x>0") {
      print(
        paste(
          sum(lower_ci > 0),
          "/",
          samples,
          " of 'x>0' 95% CIs exclude 0.",
          sep = ""
        )
      )
    }

    if (current == "x<0") {
      print(
        paste(
          sum(upper_ci < 0),
          "/",
          samples,
          " of 'x<0' 95% CIs exclude 0.",
          sep = ""
        )
      )
    }
  }
}

#' Title
#'
#' @param dataset
#' @param prefit
#'
#' @return
#'
#' @examples
parallel_run <- function(dataset, prefit, hypothesis_list, brms_formula) {
  fit <- stats::update(prefit,
    newdata = dataset,
    formula. = brms_formula,
    refresh = 0,
    silent = 2,
    chains = 4,
    control = list(adapt_delta = 0.9),
    prior = c(
      brms::prior("", class = "Intercept")
    )
  )

  result_list <- vector(mode = "list", length = length(hypothesis_list))
  names(result_list) <- hypothesis_list
  for (i in seq_along(hypothesis_list)) {
    current <- hypothesis_list[[i]]
    result_list[[current]] <- brms::hypothesis(fit, hypothesis = current)
  }
  return(result_list)
}

#' Title
#'
#' @param data_list
#' @param brms_formular
#' @param ncores
#' @param brms_family
#' @param stanvars
#' @param exports
#'
#' @return
#' @export
#'
#' @examples
calibrate_formula <- function(data_list,
                              brms_formular,
                              ncores,
                              brms_family,
                              stanvars,
                              exports,
                              hypothesis_list,
                              prefit) {
  cluster <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cluster)
  `%dopar%` <- foreach::`%dopar%`

  results <- foreach::foreach(
    dataset = data_list,
    .packages = c("brms"),
    .export = exports
  ) %dopar% {
    parallel_run(dataset, prefit, hypothesis_list, brms_formular)
  }
  parallel::stopCluster(cluster)

  result_list <- vector(mode = "list", length = length(hypothesis_list))
  names(result_list) <- hypothesis_list
  for (i in seq_along(hypothesis_list)) {
    current <- hypothesis_list[[i]]
    result_list[[current]] <- lapply(
      results,
      function(x) {
        x[[current]]
      }
    )
  }

  return(result_list)
}

#' Title
#'
#' @param dataset_N
#' @param ncores
#' @param brms_family
#' @param stanvars
#' @param exports
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
calibrate_family <- function(dataset_N,
                             ncores,
                             brms_family,
                             stanvars,
                             exports,
                             ...) {
  formula_list <- list(
    "y ~ x + z1 + z2",
    "y ~ x + z2",
    "y ~ x + z1",
    "y ~ x + z1 + z2 + z3",
    "y ~ x + z1 + z2 + z4"
  )
  data_list <- vector("list", length = dataset_N)
  for (i in seq_len(dataset_N)) {
    data_list[i] <- list(basedag_data(...))
  }
  hist(data_list[[1]]$y, main = "data", xlab = "data")

  dots <- list(...)
  dots[["x_y_coef"]] <- 0

  data_list_zero_x <- vector("list", length = dataset_N)
  for (i in seq_len(dataset_N)) {
    data_list_zero_x[i] <- list(do.call(basedag_data, dots))
  }
  hist(data_list_zero_x[[1]]$y, main = "data", xlab = "data")

  prefit <- brms::brm(
    brms::brmsformula(formula_list[[1]]),
    data = data_list[[1]],
    family = brms_family,
    stanvars = stanvars,
    chains = 0,
    refresh = 0,
    silent = 2,
    control = list(adapt_delta = 0.9),
    prior = c(
      brms::prior("", class = "Intercept")
    )
  )

  for (i in seq_along(formula_list)) {
    model_sim <- calibrate_formula(
      data_list,
      brms_formular = brms::brmsformula(formula_list[[i]]),
      ncores = ncores,
      brms_family = brms_family,
      stanvars = stanvars,
      exports = exports,
      hypothesis_list = c("x>0", "x=0"),
      prefit
    )
    calibration_result_text(
      model_name = formula_list[[i]],
      results = model_sim,
      samples = dataset_N,
      hypothesis_list = c("x>0", "x=0")
    )

    model_zero_sim <- calibrate_formula(
      data_list_zero_x,
      brms_formular = brms::brmsformula(y ~ x + z1 + z2),
      ncores = ncores,
      brms_family = brms_family,
      stanvars = stanvars,
      exports = exports,
      hypothesis_list = c("x=0", "x>0"),
      prefit
    )
    calibration_result_text(
      model_name = paste("Zero x effect"),
      results = model_zero_sim,
      samples = dataset_N,
      hypothesis_list = c("x=0", "x>0")
    )
    print("=================================================")
  }
}
