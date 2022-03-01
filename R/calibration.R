#' Title
#'
#' @param prefit
#' @param brms_formula
#' @param dataset
#' @param metric_list
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
parallel_run <- function(dataset,
                         prefit,
                         brms_formula,
                         metric_list,
                         ...) {
  metric_list <- metric_list[[1]]

  fit <- stats::update(prefit,
    newdata = dataset,
    formula. = brms::brmsformula(brms_formula),
    refresh = 0,
    silent = 2,
    chains = 4,
    inits = 1,
    control = list(adapt_delta = 0.9),
    backend = "cmdstanr",
    prior = c(
      brms::prior("", class = "Intercept")
    )
  )

  result <- metric_list_handler(fit, metric_list, ...)
  result <- append(
    result,
    list(
      "id" = dataset$id[[1]],
      "formula" = brms_formula
    )
  )
  result <- append(result, list(...))

  return(result)
}

#' Title
#'
#' @param data_list
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
calibrate_formula <- function(data_list,
                              ...) {
  `%dopar%` <- foreach::`%dopar%`

  results <- foreach::foreach(
    dataset = data_list,
    .packages = c("brms", "bayesim")
  ) %dopar% {
    parallel_run(dataset, ...)
  }
  # results <- vector(mode = "list", length = length(data_list))
  # for(i in seq_along(data_list)) {
  #   results[[i]] <- parallel_run(prefit, brms_formular, data_list[[i]], ...)
  # }

  return(do.call(rbind, lapply(results, data.frame, check.names = FALSE)))
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
                             ...) {
  data_list <- vector("list", length = dataset_N)
  for (i in seq_len(dataset_N)) {
    data <- basedag_data(...)
    data$id <- rep(i, nrow(data))
    data_list[i] <- list(data)
  }
  hist(data_list[[1]]$y, main = "data", xlab = "data")
  formula_list <- formula_list[[1]]
  combined_results_list <- vector(mode = "list", length = length(formula_list))
  for (i in seq_along(formula_list)) {
    result_df <- do.call(
      calibrate_formula,
      c(
        list(
          data_list = data_list,
          brms_formula = formula_list[[i]]
        ),
        list(...)
      )
    )
    result_df$formula <- formula_list[[i]]
    combined_results_list[[i]] <- result_df
  }
  combined_results_df <- do.call(rbind, combined_results_list)
  return(combined_results_df)
}

#' Title
#'
#' @param parameter_df
#'
#' @return
#' @export
#'
#' @examples
calibration_search <- function(parameter_df) {
  final_result <- NULL

  cluster <- parallel::makeCluster(parameter_df$ncores)
  doParallel::registerDoParallel(cluster)

  family_list <- unique(parameter_df$fit_family)

  for (f in seq_along(family_list)) {
    family <- family_list[[f]]
    family_df <- subset(parameter_df, family == family)
    link_list <- unique(family_df$fit_link)

    for (l in seq_along(link_list)) {
      link <- link_list[[l]]
      family_link_df <- subset(family_df, fit_link == link)
      prefit <- get_prefit(likelihood_lookup(
        identifier = as.character(family),
        link = as.character(link)
      ))

      for (i in seq_len(nrow(family_link_df))) {
        parameters <- as.list(family_link_df[i, ])

        row_result <- do.call(
          calibrate_family,
          c(
            parameters,
            list(
              prefit = prefit
            )
          )
        )
        final_result <- rbind(final_result, row_result)
      }
    }
  }
  parallel::stopCluster(cluster)
  return(final_result)
}
