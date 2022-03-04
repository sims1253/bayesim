#' Title
#'
#' @param dataset
#' @param metrics
#' @param seed
#' @param data_gen_conf
#' @param prefit
#' @param fit_conf
#'
#' @return
#' @export
#'
#' @examples
fit_sim <- function(prefit,
                    dataset,
                    metrics,
                    data_gen_conf,
                    fit_conf,
                    seed) {
  fit <- stats::update(prefit,
    newdata = dataset,
    formula. = brms::brmsformula(fit_conf$formula),
    refresh = 0,
    silent = 2,
    warmup = 500,
    iter = 2500,
    chains = 2,
    backend = "cmdstanr",
    seed = seed,
    init = 0.1
  )

  result <- do.call(
    metric_list_handler,
    c(
      list(fit = fit, metric_list = metrics),
      data_gen_conf
    )
  )
  result <- cbind(result, fit_conf, data_gen_conf)
  result$stan_seed <- seed
  return(as.data.frame(result))
}

#' Title
#'
#' @param metrics
#' @param seed
#' @param fit_confs
#' @param prefits
#' @param data_gen_conf
#'
#' @return
#' @export
#'
#' @examples
dataset_sim <- function(data_gen_conf,
                        fit_confs,
                        prefits,
                        metrics,
                        seed) {
  final_result <- NULL
  set.seed(seed)
  seed_list <- sample(1000000000:.Machine$integer.max,
    size = nrow(fit_confs)
  )
  datagen_result <- do.call(
    basedag_data,
    c(list(seed), data_gen_conf)
  )
  dataset <- datagen_result$dataset
  n_resample <- datagen_result$n_resample

  for (i in seq_len(nrow(fit_confs))) {
    fit_conf <- fit_confs[i, ]
    prefit <- prefits[[paste0(fit_conf$fit_family, fit_conf$fit_link)]]
    row_result <- fit_sim(
      prefit = prefit,
      dataset = dataset,
      metrics = metrics,
      data_gen_conf = data_gen_conf,
      fit_conf = fit_conf,
      seed = seed_list[[i]]
    )
    final_result <- rbind(final_result, row_result)
  }
  final_result$dataset_seed <- seed
  final_result$n_resample <- n_resample
  return(final_result)
}


#' Title
#'
#' @param data_gen_conf
#' @param seed
#' @param path
#' @param fit_confs
#' @param metrics
#' @param prefits
#'
#' @return
#' @export
#'
#' @examples
dataset_conf_sim <- function(data_gen_conf,
                             fit_confs,
                             metrics,
                             prefits,
                             seed = NULL,
                             path = NULL) {
  set.seed(seed)
  seed_list <- sample(1000000000:.Machine$integer.max,
    size = data_gen_conf$dataset_N
  )
  final_result <- NULL

  `%dopar%` <- foreach::`%dopar%`
  results <- foreach::foreach(
    seed = seed_list,
    .packages = c("brms", "bayesim")
  ) %dopar% {
    dataset_sim(
      data_gen_conf = data_gen_conf,
      fit_confs = fit_confs,
      prefits = prefits,
      metrics = metrics,
      seed = seed
    )
  }

  final_result <- do.call(rbind, results)
  final_result$data_config_seed <- seed

  if (!is.null(path)) {
    saveRDS(final_result, paste0(paste(path, data_gen_conf$id, sep = "/"), ".RDS"))
  }
  return(final_result)
}


#' Title
#'
#' @param ncores
#' @param seed
#' @param path
#' @param data_gen_confs
#' @param fit_confs
#' @param metrics
#'
#' @return
#' @export
#'
#' @examples
full_simulation <- function(data_gen_confs,
                            fit_confs,
                            metrics,
                            ncores,
                            seed = NULL,
                            path = NULL) {
  # Set seed for reproducability.
  if (!is.null(seed)) {
    set.seed(seed)
  }
  seed_list <- sample(1000000000:.Machine$integer.max,
    size = nrow(data_gen_confs)
  )
  # Multiprocessing setup TODO serial option.
  cluster <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cluster)

  # Compile a list of model configurations to be updated throughout the simulation
  # This prevents unnecessary compilation times and prevents dll overflow.
  prefit_list <- build_prefit_list(fit_configuration = fit_confs)
  final_result <- NULL

  # Iterate over dataset configurations and combine the results
  for (i in seq_len(nrow(data_gen_confs))) {
    row_result <- dataset_conf_sim(
      data_gen_conf = as.list(data_gen_confs[i, ]),
      fit_confs = fit_confs,
      metrics = metrics,
      prefits = prefit_list,
      seed = seed_list[[i]],
      path = path
    )
    final_result <- rbind(final_result, row_result)
  }
  final_result$global_seed <- seed
  # Teardown of multiprocessing setup
  parallel::stopCluster(cluster)

  if (!is.null(path)) {
    saveRDS(final_result, paste(path, "full_sim_result.RDS", sep = "/"))
  }
  return(final_result)
}


#' Title
#'
#' @param result
#'
#' @return
#' @export
#'
#' @examples
reproduce_result <- function(result) {
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
    x_y_coef = result$x_y_coef,
    y_intercept = result$y_intercept,
    sigma_y = result$sigma_y
  )

  datagen_result <- do.call(
    basedag_data,
    c(list(result$dataset_seed), data_gen_conf)
  )
  dataset <- datagen_result$dataset
  n_resample <- datagen_result$n_resample

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

  fit <- stats::update(
    prefit,
    newdata = dataset,
    formula. = brms::brmsformula(result$formula),
    refresh = 0,
    silent = 2,
    warmup = 500,
    iter = 2500,
    chains = 2,
    backend = "cmdstanr",
    seed = result$stan_seed,
    init = 0.1
  )

  return(
    list(
      fit = fit,
      dataset = dataset
    )
  )
}
