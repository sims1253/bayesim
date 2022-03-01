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
    ),
    backend = "cmdstanr"
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
  `%dopar%` <- foreach::`%dopar%`

  results <- foreach::foreach(
    dataset = data_list,
    .packages = c("brms", "bayesim")
  ) %dopar% {
    parallel_run(prefit, brms_formular, dataset, ...)
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
                             prefit,
                             ...) {
  data_list <- vector("list", length = dataset_N)
  for (i in seq_len(dataset_N)) {
    data <- basedag_data(...)
    data$id <- rep(i, nrow(data))
    data_list[i] <- list(data)
  }
  hist(data_list[[1]]$y, main = "data", xlab = "data")

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

#' Title
#'
#' @param result_df
#' @param formula_list
#'
#' @return
#' @export
#'
#' @examples
viz_calibration <- function(result_df, formula_list, title, y_lab) {
  formula_list <- forcats::fct_inorder(unlist(simulation_parameters$formula_list))
  eq_list <- vector(mode = "numeric", length = length(formula_list))
  greater_list <- vector(mode = "numeric", length = length(formula_list))
  for (i in seq_along(formula_list)) {
    current <- formula_list[[i]]
    subset_df <- dplyr::filter(result_df, formula == current)
    greater_list[[i]] <- sum(subset_df[["hyp_x>0_lower_ci"]] > 0) / nrow(subset_df)
    eq_list[[i]] <- sum(subset_df[["hyp_x=0_lower_ci"]] > 0 | subset_df[["hyp_x=0_upper_ci"]] < 0) / nrow(subset_df)
  }
  df <- data.frame(formula = formula_list, "x=0" = eq_list, "x>0" = greater_list, check.names = FALSE)

  long <- reshape2::melt(df, c("formula"))

  ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = formula,
      y = value,
      colour = formula,
      shape = formula
    )
  ) +
    ggplot2::geom_point(size = 5) +
    ggplot2::facet_wrap("variable") +
    ggplot2::scale_colour_viridis_d() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::ggtitle(title) +
    ggplot2::ylab(y_lab)
}

#' Title
#'
#' @param x_y_coef_list
#' @param y_intercept_list
#' @param link_list
#' @param likelihood_list
#' @param sigma_y_list
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
calibration_search <- function(x_y_coef_list,
                               y_intercept_list,
                               link_list,
                               likelihood_list,
                               sigma_y_list,
                               search_list,
                               ncores,
                               ...) {
  final_result <- NULL
  cluster <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cluster)

  for (l in seq_along(likelihood_list)) {
    likelihood <- likelihood_list[[l]]

    for (k in seq_along(link_list)) {
      link <- link_list[[k]]
      prefit <- get_prefit(likelihood_lookup(
        identifier = likelihood,
        link = link
      ))
      y_intercept_list <- search_list[[likelihood]][[link]]$y_intercept_list
      x_y_coef_list <- search_list[[likelihood]][[link]]$x_y_coef_list

      for (i in seq_along(x_y_coef_list)) {
        for (j in seq_along(y_intercept_list)) {
          if (is.null(final_result)) {
            result_df <- calibrate_family(
              x_y_coef = x_y_coef_list[[i]],
              y_intercept = y_intercept_list[[j]],
              data_link = brms:::inv_link(link),
              RNG = rng_lookup(likelihood),
              sigma_y = sigma_y_list[[likelihood]],
              prefit = prefit,
              ...
            )
            result_df$link <- link
            result_df$family <- likelihood

            final_result <- result_df
          } else {
            result_df <- calibrate_family(
              x_y_coef = x_y_coef_list[[i]],
              y_intercept = y_intercept_list[[j]],
              data_link = brms:::inv_link(link),
              RNG = rng_lookup(likelihood),
              sigma_y = sigma_y_list[[likelihood]],
              prefit = prefit,
              ...
            )
            result_df$link <- link
            result_df$family <- likelihood

            final_result <- rbind(final_result, result_df)
          }
        }
      }
    }
  }
  parallel::stopCluster(cluster)
  return(final_result)
}
