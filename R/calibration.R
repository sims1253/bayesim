#' Title
#'
#' @param result_df
#' @param formula_list
#'
#' @return
#' @export
#'
#' @examples
calibration_text_result <- function(title,
                                    result_df,
                                    formula_list) {
  print(title)
  for (i in seq_along(formula_list)) {
    current <- formula_list[[i]]
    subset_df <- dplyr::filter(result_df, formula == current)

    print(current)
    print(
      paste0(
        "The mean posterior bias was ",
        mean(subset_df$bias),
        "."
      )
    )

    print(
      paste0(
        "The mean posterior rmse was ",
        mean(subset_df$rmse_s),
        "."
      )
    )

    print(
      paste0(
        "The mean confidence interval width was ",
        mean(subset_df[["hyp_x>0_upper_ci"]] - subset_df[["hyp_x>0_lower_ci"]]),
        "."
      )
    )

    print(
      paste0(
        sum(subset_df[["hyp_x>0_lower_ci"]] > 0),
        "/",
        nrow(subset_df),
        " x>0 90% CIs exclude 0."
      )
    )

    print(
      paste0(
        sum(subset_df[["hyp_x=0_lower_ci"]] > 0 | subset_df[["hyp_x=0_upper_ci"]] < 0),
        "/",
        nrow(subset_df),
        " x=0 95% CIs exclude 0."
      )
    )
    print("====================================")
  }
}


#' Title
#'
#' @param prefit
#' @param brms_formula
#' @param dataset
#' @param metric_list
#' @param ...
#'
#' @return
#'
#' @examples
parallel_run <- function(prefit,
                         brms_formula,
                         dataset,
                         metric_list,
                         data_link,
                         RNG,
                         ...) {
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

  result <- metric_list_handler(fit, metric_list, ...)
  result <- append(result, list("id" = dataset$id[[1]]))
  result <- append(result, list(...))

  return(result)
}


#' Title
#'
#' @param data_list
#' @param brms_formular
#' @param prefit
#' @param ncores
#' @param ...
#'
#' @return
#'
#' @examples
calibrate_formula <- function(data_list,
                              brms_formular,
                              prefit,
                              ncores,
                              ...) {
  cluster <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cluster)
  `%dopar%` <- foreach::`%dopar%`

  results <- foreach::foreach(
    dataset = data_list,
    .packages = c("brms", "bayesim")
  ) %dopar% {
    parallel_run(prefit, brms_formular, dataset, ...)
  }
  parallel::stopCluster(cluster)
  # results <- vector(mode = "list", length = length(data_list))
  # for(i in seq_along(data_list)) {
  #   results[[i]] <- parallel_run(prefit, brms_formular, data_list[[i]], ...)
  # }

  return(do.call(rbind,lapply(results,data.frame, check.names = FALSE)))
}


#' Title
#'
#' @param dataset_N
#' @param formula_list
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
calibrate_family <- function(dataset_N,
                             formula_list,
                             brms_family = brms_family,
                             stanvars = stanvars,
                             ...) {
  data_list <- vector("list", length = dataset_N)
  for (i in seq_len(dataset_N)) {
    data <- basedag_data(...)
    data$id <- rep(i, nrow(data))
    data_list[i] <- list(data)
  }
  hist(data_list[[1]]$y, main = "data", xlab = "data")

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

  combined_results_list <- vector(mode = "list", length = length(formula_list))
  for (i in seq_along(formula_list)) {
    result_df <- calibrate_formula(
      data_list,
      brms_formular = brms::brmsformula(formula_list[[i]]),
      prefit = prefit,
      ...
    )
    result_df$formula <- rep(formula_list[[i]], nrow(result_df))
    combined_results_list[[i]] <- result_df
  }
  combined_results_df <- do.call(rbind, combined_results_list)
}
