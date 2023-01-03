#' Title
#'
#' @param dataset
#' @param seed
#' @param data_gen_conf
#' @param prefit
#' @param fit_conf
#' @param testing_data
#' @param brms_backend
#' @param debug
#' @param path
#'
#' @return
#' @export
#'
#' @examples
fit_sim <- function(prefit,
                    dataset,
                    data_gen_output,
                    fit_conf,
                    seed,
                    debug,
                    result_path,
                    stan_pars,
                    ...) {
  fit <- stats::update(prefit,
    newdata = dataset,
    formula. = brms::brmsformula(fit_conf$formula),
    refresh = 0,
    silent = 2,
    warmup = stan_pars$warmup,
    iter = stan_pars$iter,
    chains = stan_pars$chains,
    backend = stan_pars$backend,
    seed = seed,
    init = stan_pars$init
  )
  if (debug == TRUE) {
    saveRDS(fit, paste0(paste(result_path, "fit", sep = "/"), ".RDS"))
  }

  all_metric_results <- do.call(
    metric_list_handler,
    c(
      list(
        fit = fit,
        data_gen_output = data_gen_output,
        fit_conf = fit_conf
      ),
      list(...)
    )
  )
  if (debug == TRUE) {
    saveRDS(all_metric_results, paste0(paste(result_path, "metric_results", sep = "/"), ".RDS"))
  }
  return(all_metric_results)
}

#' Title
#'
#' @param seed
#' @param fit_confs
#' @param prefits
#' @param data_gen_conf
#' @param brms_backend
#' @param cmdstan_path
#' @param seed
#' @param debug
#' @param path
#'
#' @return
#' @export
#'
#' @examples
dataset_sim <- function(data_gen_conf,
                        data_gen_fun,
                        fit_confs,
                        prefits,
                        stan_pars,
                        seed,
                        debug,
                        result_path,
                        ...) {
  if (stan_pars$backend == "cmdstanr") {
    cmdstanr::set_cmdstan_path(stan_pars$cmdstan_path)
    if (!is.null(stan_pars$cmdstan_write_path)) {
      options(cmdstanr_write_stan_file_dir = stan_pars$cmdstan_write_path)
    }
  }
  final_result <- vector(mode = "list", length = nrow(fit_confs))
  loo_objects <- vector(mode = "list", length = nrow(fit_confs))
  set.seed(seed)
  seed_list <- sample(1000000000:.Machine$integer.max,
    size = nrow(fit_confs)
  )
  datagen_result <- do.call(
    data_gen_fun,
    c(list(seed = seed), data_gen_conf)
  )
  if (debug == TRUE) {
    saveRDS(datagen_result, paste0(paste(result_path, "datagen_result", sep = "/"), ".RDS"))
    saveRDS(data_gen_conf, paste0(paste(result_path, "data_gen_conf", sep = "/"), ".RDS"))
  }

  for (i in seq_len(nrow(fit_confs))) {
    fit_conf <- fit_confs[i, ]
    prefit <- prefits[[fit_conf_key(fit_conf)]]
    if (debug == TRUE) {
      saveRDS(fit_conf, paste0(paste(result_path, "fit_conf", sep = "/"), ".RDS"))
      saveRDS(prefit, paste0(paste(result_path, "prefit", sep = "/"), ".RDS"))
    }

    final_result[[i]] <- fit_sim(
      prefit = prefit,
      dataset = datagen_result$dataset,
      testing_data = datagen_result$testing_data,
      data_gen_output = datagen_result$data_gen_output,
      fit_conf = fit_conf,
      seed = seed_list[[i]],
      debug = debug,
      result_path = result_path,
      stan_pars = stan_pars,
      ...
    )
  }

  final_result <- dplyr::bind_rows(final_result)
  if ("NA." %in% colnames(final_result)) {
    final_result <- subset(final_result,
      select = -c(which(colnames(final_result) == "NA."))
    )
  }

  final_result$dataset_seed <- seed

  if (debug == TRUE) {
    saveRDS(final_result, paste0(paste(result_path, "dataset_result", sep = "/"), ".RDS"))
  }

  return(dplyr::as_tibble(final_result))
}


#' Title
#'
#' @param data_gen_conf
#' @param seed
#' @param path
#' @param fit_confs
#' @param prefits
#' @param brms_backend
#' @param ncores
#' @param debug
#'
#' @return
#' @export
#'
#' @examples
dataset_conf_sim <- function(data_gen_conf,
                             data_gen_fun,
                             fit_confs,
                             prefits,
                             seed = NULL,
                             result_path = NULL,
                             stan_pars,
                             ncores,
                             cluster_type,
                             debug,
                             global_seed,
                             ...) {
  set.seed(seed)
  seed_list <- sample(1000000000:.Machine$integer.max,
    size = data_gen_conf$dataset_N
  )

  if (file.exists(paste0(paste(result_path, data_gen_conf$id, sep = "/"), ".RDS"))) {
    return(readRDS(paste0(paste(result_path, data_gen_conf$id, sep = "/"), ".RDS")))
  } else {
    if (ncores > 1) {
      # Multiprocessing setup
      cluster <- parallel::makeCluster(ncores,
        type = cluster_type,
        outfile = paste(
          result_path, "cluster_log",
          sep = "/"
        )
      )
      doParallel::registerDoParallel(cluster)
      parallel::clusterEvalQ(cl = cluster, {
        options(mc.cores = 1)
      })
      `%dopar%` <- foreach::`%dopar%`

      # Multiprocessing run
      final_result <- foreach::foreach(
        par_seed = seed_list,
        .packages = "bayesim",
        .combine = dplyr::bind_rows, .multicombine = TRUE,
        .maxcombine = length(seed_list),
        .verbose = debug,
        .inorder = FALSE
      ) %dopar% {
        dataset_sim(
          data_gen_conf = data_gen_conf,
          data_gen_fun = data_gen_fun,
          fit_confs = fit_confs,
          prefits = prefits,
          stan_pars = stan_pars,
          seed = par_seed,
          debug = debug,
          result_path = result_path,
          ...
        )
      }

      # Multiprocessing teardown
      parallel::stopCluster(cluster)
    } else {
      results <- vector(mode = "list", length = length(seed_list))
      for (i in seq_along(seed_list)) {
        results[[i]] <- dataset_sim(
          data_gen_conf = data_gen_conf,
          data_gen_fun = data_gen_fun,
          fit_confs = fit_confs,
          prefits = prefits,
          stan_pars = stan_pars,
          seed = seed_list[[i]],
          debug = debug,
          result_path = result_path,
          ...
        )
      }
      final_result <- dplyr::bind_rows(results)
    }

    final_result$data_config_seed <- seed
    final_result$global_seed <- global_seed
    final_result$brms_backend <- stan_pars$backend

    if (!is.null(result_path)) {
      saveRDS(
        final_result,
        paste0(paste(result_path, data_gen_conf$id, sep = "/"), ".RDS")
      )
    }
    return(final_result)
  }
}


#' Title
#'
#' @param seed
#' @param path
#' @param data_gen_confs
#' @param ncores_prefit
#' @param ncores_simulation
#' @param brms_backend
#' @param fit_confs
#' @param debug
#'
#' @return
#' @export
#'
#' @examples
full_simulation <- function(data_gen_confs,
                            data_gen_fun = NULL,
                            fit_confs,
                            ncores_simulation = 1,
                            cluster_type = "PSOCK",
                            stan_pars,
                            seed = NULL,
                            result_path = NULL,
                            debug = FALSE,
                            ...) {
  # Set seed for reproducability.
  if (!is.null(seed)) {
    set.seed(seed)
  }
  seed_list <- sample(1000000000:.Machine$integer.max,
    size = nrow(data_gen_confs)
  )
  if (stan_pars$backend == "cmdstanr") {
    cmdstanr::set_cmdstan_path(stan_pars$cmdstan_path)
    if (!is.null(stan_pars$cmdstan_write_path)) {
      options(cmdstanr_write_stan_file_dir = stan_pars$cmdstan_write_path)
    }
  }

  # Compile a list of model configurations to be updated throughout the run
  # This prevents unnecessary compilation times and prevents dll overflow.
  prefit_list <- build_prefit_list(
    fit_configuration = fit_confs,
    stan_pars = stan_pars
  )
  final_result <- vector(mode = "list", length = nrow(data_gen_confs))

  # Iterate over dataset configurations and combine the results
  for (i in seq_len(nrow(data_gen_confs))) {
    final_result[[i]] <- dataset_conf_sim(
      data_gen_conf = as.list(data_gen_confs[i, ]),
      data_gen_fun = data_gen_fun,
      fit_confs = fit_confs,
      prefits = prefit_list,
      seed = seed_list[[i]],
      result_path = result_path,
      stan_pars = stan_pars,
      ncores = ncores_simulation,
      cluster_type = cluster_type,
      debug = debug,
      global_seed = seed,
      ...
    )
  }
  final_result <- dplyr::bind_rows(final_result)

  if (!is.null(result_path)) {
    saveRDS(final_result, paste(result_path, "full_sim_result.RDS", sep = "/"))
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
reproduce_result <- function(result, data_gen_fun) {
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
    backend = result$brms_backend,
    prior = prior_lookup(result$fit_family)
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
    shape = result$shape,
    seed = result$dataset_seed
  )

  datagen_result <- do.call(
    data_gen_fun,
    data_gen_conf
  )
  dataset <- datagen_result$dataset
  sampling_loops <- datagen_result$sampling_loops
  bad_samples <- datagen_result$bad_samples
  testing_data <- datagen_result$testing_data

  fit <- stats::update(prefit,
    newdata = dataset,
    formula. = brms::brmsformula(result$formula),
    refresh = 0,
    silent = 2,
    warmup = 500,
    iter = 2500,
    chains = 2,
    backend = result$brms_backend,
    seed = result$stan_seed,
    init = 0.1
  )

  return(
    list(
      fit = fit,
      dataset = dataset,
      testing_data = datagen_result$testing_data,
      sampling_loops = sampling_loops,
      bad_samples = bad_samples,
      data_gen_conf = data_gen_conf
    )
  )
}
