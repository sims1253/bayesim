#' Title
#'
#' @param dataset
#' @param seed
#' @param data_gen_conf
#' @param prefit
#' @param fit_conf
#' @param testing_data
#' @param numeric_metrics
#' @param predictive_metrics
#' @param brms_backend
#'
#' @return
#' @export
#'
#' @examples
fit_sim <- function(prefit,
                    dataset,
                    testing_data,
                    numeric_metrics,
                    predictive_metrics,
                    data_gen_conf,
                    fit_conf,
                    brms_backend,
                    seed) {
  fit <- stats::update(prefit,
    newdata = dataset,
    formula. = brms::brmsformula(fit_conf$formula),
    refresh = 0,
    silent = 2,
    warmup = 500,
    iter = 2500,
    chains = 2,
    backend = brms_backend,
    seed = seed,
    init = 0.1
  )

  all_metric_results <- do.call(
    metric_list_handler,
    c(
      list(
        fit = fit,
        numeric_metrics = numeric_metrics,
        predictive_metrics = predictive_metrics,
        testing_data = testing_data
      ),
      data_gen_conf
    )
  )
  numeric_results <- all_metric_results$numeric_results
  loo_objects <- all_metric_results$loo_objects
  final_result <- list(
    numeric_results = data.frame(
      c(
        numeric_results,
        fit_conf,
        data_gen_conf,
        c(stan_seed = seed)
      )
    ),
    loo_objects = loo_objects
  )
  return(final_result)
}

#' Title
#'
#' @param seed
#' @param fit_confs
#' @param prefits
#' @param data_gen_conf
#' @param numeric_metrics
#' @param predictive_metrics
#' @param brms_backend
#'
#' @return
#' @export
#'
#' @examples
dataset_sim <- function(data_gen_conf,
                        fit_confs,
                        prefits,
                        numeric_metrics,
                        predictive_metrics,
                        brms_backend,
                        seed) {
  final_result <- vector(mode = "list", length = nrow(fit_confs))
  loo_objects <- vector(mode = "list", length = nrow(fit_confs))
  set.seed(seed)
  seed_list <- sample(1000000000:.Machine$integer.max,
    size = nrow(fit_confs)
  )
  datagen_result <- do.call(
    basedag_data,
    c(list(seed = seed), data_gen_conf)
  )
  dataset <- datagen_result$dataset
  sampling_loops <- datagen_result$sampling_loops
  bad_samples <- datagen_result$bad_samples
  testing_data <- datagen_result$testing_data

  for (i in seq_len(nrow(fit_confs))) {
    fit_conf <- fit_confs[i, ]
    prefit <- prefits[[paste0(fit_conf$fit_family, fit_conf$fit_link)]]
    row_results <- fit_sim(
      prefit = prefit,
      dataset = dataset,
      testing_data = testing_data,
      numeric_metrics = numeric_metrics,
      predictive_metrics = predictive_metrics,
      data_gen_conf = data_gen_conf,
      fit_conf = fit_conf,
      brms_backend = brms_backend,
      seed = seed_list[[i]]
    )
    final_result[[i]] <- row_results$numeric_results
    loo_objects[[i]] <- row_results$loo_objects
  }

  names(loo_objects) <- seq_len(length(loo_objects))
  loo_compare_results <- loo_compare_handler(loo_objects, predictive_metrics)

  final_result <- do.call(rbind, final_result) # TODO handle predictive performance results
  final_result <- cbind(final_result, loo_compare_results)

  final_result$dataset_seed <- seed
  final_result$bad_samples <- bad_samples
  final_result$sampling_loops <- sampling_loops
  return(as.data.frame(final_result))
}


#' Title
#'
#' @param data_gen_conf
#' @param seed
#' @param path
#' @param fit_confs
#' @param prefits
#' @param numeric_metrics
#' @param predictive_metrics
#' @param brms_backend
#' @param ncores
#'
#' @return
#' @export
#'
#' @examples
dataset_conf_sim <- function(data_gen_conf,
                             fit_confs,
                             numeric_metrics,
                             predictive_metrics,
                             prefits,
                             seed = NULL,
                             path = NULL,
                             brms_backend,
                             ncores) {
  set.seed(seed)
  seed_list <- sample(1000000000:.Machine$integer.max,
    size = data_gen_conf$dataset_N
  )

  if (file.exists(paste0(paste(path, data_gen_conf$id, sep = "/"), ".RDS"))) {
    return(readRDS(paste0(paste(path, data_gen_conf$id, sep = "/"), ".RDS")))
  } else {
    if (ncores > 1) {
      # Multiprocessing setup
      cluster <- parallel::makeCluster(ncores, type = "PSOCK")
      doParallel::registerDoParallel(cluster)
      parallel::clusterEvalQ(cl = cluster, {
        library(brms)
        library(bayesim)
        options(mc.cores = 1)
      })
      `%dopar%` <- foreach::`%dopar%`

      # Multiprocessing run
      results <- foreach::foreach(
        par_seed = seed_list
      ) %dopar% {
        dataset_sim(
          data_gen_conf = data_gen_conf,
          fit_confs = fit_confs,
          prefits = prefits,
          numeric_metrics = numeric_metrics,
          predictive_metrics = predictive_metrics,
          brms_backend = brms_backend,
          seed = par_seed
        )
      }

      # Multiprocessing teardown
      parallel::stopCluster(cluster)
    } else {
      results <- vector(mode = "list", length = length(seed_list))
      for (i in seq_along(seed_list)) {
        results[[i]] <- dataset_sim(
          data_gen_conf = data_gen_conf,
          fit_confs = fit_confs,
          prefits = prefits,
          numeric_metrics = numeric_metrics,
          predictive_metrics = predictive_metrics,
          brms_backend = brms_backend,
          seed = seed_list[[i]]
        )
      }
    }

    final_result <- do.call(rbind, results)
    final_result$data_config_seed <- seed

    if (!is.null(path)) {
      saveRDS(final_result, paste0(paste(path, data_gen_conf$id, sep = "/"), ".RDS"))
    }
    return(final_result)
  }
}


#' Title
#'
#' @param seed
#' @param path
#' @param data_gen_confs
#' @param numeric_metrics
#' @param predictive_metrics
#' @param ncores_prefit
#' @param ncores_simulation
#' @param brms_backend
#' @param fit_confs
#'
#' @return
#' @export
#'
#' @examples
full_simulation <- function(data_gen_confs,
                            fit_confs,
                            numeric_metrics,
                            predictive_metrics,
                            ncores_simulation = 1,
                            brms_backend = "cmdstanr",
                            seed = NULL,
                            path = NULL) {
  # Set seed for reproducability.
  if (!is.null(seed)) {
    set.seed(seed)
  }
  seed_list <- sample(1000000000:.Machine$integer.max,
    size = nrow(data_gen_confs)
  )

  # Compile a list of model configurations to be updated throughout the simulation
  # This prevents unnecessary compilation times and prevents dll overflow.
  prefit_list <- build_prefit_list(fit_configuration = fit_confs)
  final_result <- vector(mode = "list", length = nrow(data_gen_confs))

  # Iterate over dataset configurations and combine the results
  for (i in seq_len(nrow(data_gen_confs))) {
    final_result[[i]] <- dataset_conf_sim(
      data_gen_conf = as.list(data_gen_confs[i, ]),
      fit_confs = fit_confs,
      numeric_metrics,
      predictive_metrics,
      prefits = prefit_list,
      seed = seed_list[[i]],
      path = path,
      brms_backend = brms_backend,
      ncores = ncores_simulation
    )
  }
  final_result <- do.call(rbind, final_result)
  final_result$global_seed <- seed
  final_result$brms_backend <- brms_backend

  if (!is.null(path)) {
    saveRDS(final_result, paste(path, "full_sim_result.RDS", sep = "/"))
  }
  return(final_result)
}


#' This method will reproduce the exact dataset and fit corresponding to the
#' supplied result dataframe row.
#'
#' The code in this function is written so that all seeds are set at the right
#' time and all following code after setting the seed replicates exactly as
#' during the simulation.
#'
#' @param result
#'
#' @return
#' @export
#'
#' @examples
reproduce_result <- function(result) {
  family <- brms_family_lookup(
    result$fit_family,
    result$fit_link
  )
  prefit <- brms::brm(
    y ~ 1 + x,
    data = list(y = c(0.5), x = c(1)),
    family = family,
    stanvars = family$stanvars,
    chains = 0,
    refresh = 0,
    silent = 2,
    backend = "cmdstanr",
    prior = c(
      brms::set_prior("", class = "Intercept"),
      brms::set_prior("", class = second_family_parameter_lookup(result$fit_family))
    ),
    init = 0.1
  )

  data_gen_conf <- list(
    z1_x_coef = result$z1_x_coef,
    z3_x_coef = result$z3_x_coef,
    z1_y_coef = result$z1_y_coef,
    z2_y_coef = result$z2_y_coef,
    x_z4_coef = result$x_z4_coef,
    y_z4_coef = result$y_z4_coef,
    sigma_z1 = result$sigma_z1,
    sigma_z2 = result$sigma_z2,
    sigma_z3 = result$sigma_z3,
    sigma_z4 = result$sigma_z4,
    sigma_x = result$sigma_x,
    data_N = result$data_N,
    dataset_N = result$dataset_N,
    data_family = result$data_family,
    data_link = result$data_link,
    lb = result$lb,
    ub = result$ub,
    resample = result$resample,
    x_y_coef = result$x_y_coef,
    y_intercept = result$y_intercept,
    sigma_y = result$sigma_y,
    shape = result$shape
  )

  datagen_result <- do.call(
    basedag_data,
    c(list(seed = seed), data_gen_conf)
  )
  dataset <- datagen_result$dataset
  sampling_loops <- datagen_result$sampling_loops
  bad_samples <- datagen_result$bad_samples
  testing_data <- datagen_result$testing_data

  fit <- stats::update(prefit,
    newdata = dataset,
    formula. = brms::brmsformula(fit_conf$formula),
    refresh = 0,
    silent = 2,
    warmup = 500,
    iter = 2500,
    chains = 2,
    backend = result$brms_backend,
    seed = seed,
    init = 0.1
  )

  return(
    list(
      fit = fit,
      dataset = dataset,
      testing_data = datagen_result$testing_data,
      sampling_loops = sampling_loops,
      bad_samples = bad_samples
    )
  )
}
