#' Fit a single model in a simulation study (DEPRECATED)
#'
#' @description
#' `r lifecycle::deprecated("fit_sim()")`
#' This function is deprecated. New workflows should use `simulation_config()` + `run_simulation()`.
#'
#' @param prefit A precompiled brms model object from \code{\link{get_prefit}}
#' @param dataset A data.frame containing the data to fit the model to
#' @param data_gen_output Output from the data generation function containing true parameters
#' @param fit_conf A list or data.frame row containing model configuration (formula, family, etc.)
#' @param seed Random seed for reproducible model fitting
#' @param debug Logical; if TRUE, saves intermediate results for debugging
#' @param result_path Path where debug files should be saved
#' @param stan_pars List of Stan/brms parameters (chains, iter, warmup, etc.)
#' @param ... Additional arguments passed to metric calculation functions
#'
#' @return A tibble containing calculated metrics for the fitted model
#' @keywords internal
#' @export
#'
#' @examples
#' \dontrun{
#' # This function is deprecated. Use simulation_config() + run_simulation() instead.
#' }
fit_sim <- function(
  prefit,
  dataset,
  data_gen_output,
  fit_conf,
  seed,
  debug,
  result_path,
  stan_pars,
  ...
) {
  fit <- stats::update(
    prefit,
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
  if (debug) {
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
  if (debug) {
    saveRDS(
      all_metric_results,
      paste0(paste(result_path, "metric_results", sep = "/"), ".RDS")
    )
  }
  all_metric_results
}

#' Run simulation for a single dataset configuration (DEPRECATED)
#'
#' @description
#' `r lifecycle::deprecated("dataset_sim()")`
#' This function is deprecated. New workflows should use `simulation_config()` + `run_simulation()`.
#'
#' @param data_gen_conf A list containing data generation configuration parameters
#' @param fit_confs A data.frame containing model fitting configurations
#' @param prefits A list of precompiled brms model objects
#' @param stan_pars A list containing Stan/brms fitting parameters
#' @param seed Random seed for reproducible data generation and model fitting
#' @param debug Logical; if TRUE, saves intermediate results for debugging
#' @param result_path Path where debug files should be saved
#' @param data_generator_fn A function with signature `(config, seed) -> list`
#'   that generates data. The returned list must contain: `dataset`,
#'   `sampling_loops`, `bad_samples`, `testing_data`, and `true_parameters`.
#' @param ... Additional arguments passed to metric calculation functions
#'
#' @return A tibble containing simulation results for all models
#' @keywords internal
#' @export
#'
#' @examples
#' \dontrun{
#' # This function is deprecated. Use simulation_config() + run_simulation() instead.
#' }
dataset_sim <- function(
  data_gen_conf,
  fit_confs,
  prefits,
  stan_pars,
  seed,
  debug,
  result_path,
  data_generator_fn = NULL,
  ...
) {
  if (stan_pars$backend == "cmdstanr") {
    cmdstanr::set_cmdstan_path(stan_pars$cmdstan_path)
    if (!is.null(stan_pars$cmdstan_write_path)) {
      options(cmdstanr_write_stan_file_dir = stan_pars$cmdstan_write_path)
    }
  }
  final_result <- vector(mode = "list", length = nrow(fit_confs))
  loo_objects <- vector(mode = "list", length = nrow(fit_confs))
  set.seed(seed)
  seed_list <- sample(1000000000:.Machine$integer.max, size = nrow(fit_confs))

  if (is.null(data_generator_fn)) {
    stop(
      "data_generator_fn must be provided. ",
      "Pass a function with signature (config, seed) -> list"
    )
  }

  datagen_result <- data_generator_fn(config = data_gen_conf, seed = seed)
  if (debug) {
    saveRDS(
      datagen_result,
      paste0(paste(result_path, "datagen_result", sep = "/"), ".RDS")
    )
    saveRDS(
      data_gen_conf,
      paste0(paste(result_path, "data_gen_conf", sep = "/"), ".RDS")
    )
  }

  for (i in seq_len(nrow(fit_confs))) {
    fit_conf <- fit_confs[i, ]
    prefit <- prefits[[fit_conf_key(fit_conf)]]
    if (debug) {
      saveRDS(
        fit_conf,
        paste0(paste(result_path, "fit_conf", sep = "/"), ".RDS")
      )
      saveRDS(prefit, paste0(paste(result_path, "prefit", sep = "/"), ".RDS"))
    }

    final_result[[i]] <- fit_sim(
      prefit = prefit,
      dataset = datagen_result$dataset,
      testing_data = datagen_result$testing_data,
      data_gen_output = datagen_result$true_parameters,
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
    final_result <- subset(
      final_result,
      select = -c(which(colnames(final_result) == "NA."))
    )
  }

  final_result$dataset_seed <- seed

  if (debug) {
    saveRDS(
      final_result,
      paste0(paste(result_path, "dataset_result", sep = "/"), ".RDS")
    )
  }

  dplyr::as_tibble(final_result)
}


#' Run simulation for a single data generation configuration (DEPRECATED)
#'
#' @description
#' `r lifecycle::deprecated("dataset_conf_sim()")`
#' This function is deprecated. New workflows should use `simulation_config()` + `run_simulation()`.
#'
#' @param data_gen_conf A list containing data generation configuration parameters
#' @param fit_confs A data.frame containing model fitting configurations
#' @param prefits A list of precompiled brms model objects
#' @param seed Random seed for reproducible data generation
#' @param result_path Path where results and debug files should be saved
#' @param stan_pars A list containing Stan/brms fitting parameters (backend, chains, iter, warmup, etc.)
#' @param ncores Number of cores to use for parallel processing
#' @param cluster_type Type of cluster for parallel processing (e.g., "PSOCK")
#' @param debug Logical; if TRUE, saves intermediate results for debugging
#' @param global_seed Global seed for the entire simulation run
#' @param ... Additional arguments passed to dataset_sim and metric calculation functions
#'
#' @return A tibble containing simulation results for all datasets and models
#' @keywords internal
#' @export
#'
#' @examples
#' \dontrun{
#' # This function is deprecated. Use simulation_config() + run_simulation() instead.
#' }
dataset_conf_sim <- function(
  data_gen_conf,
  fit_confs,
  prefits,
  seed = NULL,
  result_path = NULL,
  stan_pars,
  ncores,
  cluster_type,
  debug,
  global_seed,
  ...
) {
  set.seed(seed)
  seed_list <- sample(
    1000000000:.Machine$integer.max,
    size = data_gen_conf$dataset_N
  )

  if (
    file.exists(
      paste0(paste(result_path, data_gen_conf$id, sep = "/"), ".RDS")
    )
  ) {
    readRDS(paste0(paste(result_path, data_gen_conf$id, sep = "/"), ".RDS"))
  } else {
    if (ncores > 1) {
      if (debug) {
        cluster <- parallel::makeCluster(
          ncores,
          type = cluster_type,
          outfile = paste(result_path, "cluster_log")
        )
      } else {
        cluster <- parallel::makeCluster(ncores, type = cluster_type)
      }
      # Multiprocessing setup
      doParallel::registerDoParallel(cluster)
      parallel::clusterEvalQ(cl = cluster, {
        options(mc.cores = 1)
      })
      `%dopar%` <- foreach::`%dopar%`

      # Multiprocessing run
      final_result <- foreach::foreach(
        par_seed = seed_list,
        .packages = "bayesim",
        .combine = dplyr::bind_rows,
        .multicombine = TRUE,
        .maxcombine = length(seed_list),
        .verbose = debug,
        .inorder = FALSE
      ) %dopar%
        {
          dataset_sim(
            data_gen_conf = data_gen_conf,
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
    final_result
  }
}


#' Run a full simulation study (DEPRECATED)
#'
#' @description
#' `r lifecycle::deprecated("full_simulation()")`
#' This function is deprecated. New workflows should use `simulation_config()` + `run_simulation()`.
#'
#' Legacy interface for older bayesim simulation studies.
#'
#' New code should prefer `simulation_config()`, `run_simulation()`, and
#' `resume_simulation()`. This wrapper is retained for compatibility with the
#' pre-rewrite API.
#'
#' @param data_gen_confs A data.frame where each row specifies a data generation configuration
#' @param fit_confs A data.frame where each row specifies a model fitting configuration
#' @param ncores_simulation Number of cores to use for parallel dataset simulation
#' @param cluster_type Type of cluster for parallel processing (e.g., "PSOCK", "FORK")
#' @param stan_pars A list containing Stan/brms fitting parameters:
#'   \itemize{
#'     \item backend: Stan backend ("cmdstanr" or "rstan")
#'     \item chains: Number of MCMC chains
#'     \item iter: Number of iterations per chain
#'     \item warmup: Number of warmup iterations
#'     \item init: Initial values for MCMC
#'     \item cmdstan_path: Path to CmdStan (if using cmdstanr backend)
#'     \item cmdstan_write_path: Path for CmdStan output files
#'   }
#' @param seed Random seed for reproducible simulation
#' @param result_path Path where results should be saved (as RDS files)
#' @param debug Logical; if TRUE, enables verbose output and saves intermediate results
#' @param ... Additional arguments passed to dataset_conf_sim and metric calculation functions
#'
#' @return A tibble containing all simulation results
#' @keywords internal
#' @export
#'
#' @examples
#' \dontrun{
#' # This function is deprecated. Use simulation_config() + run_simulation() instead.
#' # Define data generation configurations
#' data_gen_confs <- data.frame(
#'   id = c("config1", "config2"),
#'   data_N = c(100, 100),
#'   dataset_N = c(10, 10)
#' )
#'
#' # Define model fitting configurations
#' fit_confs <- data.frame(
#'   fit_family = c("gaussian", "student_t"),
#'   fit_link = c("identity", "identity"),
#'   formula = c("y ~ x", "y ~ x")
#' )
#'
#' # Set Stan parameters
#' stan_pars <- list(
#'   backend = "cmdstanr",
#'   chains = 4,
#'   iter = 2000,
#'   warmup = 1000,
#'   init = 0.1
#' )
#'
#' # Run simulation
#' results <- full_simulation(
#'   data_gen_confs = data_gen_confs,
#'   fit_confs = fit_confs,
#'   stan_pars = stan_pars,
#'   ncores_simulation = 4,
#'   seed = 12345
#' )
#' }
full_simulation <- function(
  data_gen_confs,
  fit_confs,
  ncores_simulation = 1,
  cluster_type = "PSOCK",
  stan_pars,
  seed = NULL,
  result_path = NULL,
  debug = FALSE,
  ...
) {
  # Set seed for reproducability.
  if (!is.null(seed)) {
    set.seed(seed)
  }
  seed_list <- sample(
    1000000000:.Machine$integer.max,
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
  final_result
}


#' This method will reproduce the exact dataset and fit corresponding to the
#' supplied result dataframe row.
#'
#' The code in this function is written so that all seeds are set at the right
#' time and all following code after setting the seed replicates exactly as
#' during the simulation.
#'
#' @param result A data.frame row containing simulation result metadata
#' @param data_generator_fn A function with signature `(config, seed) -> list`
#'   that generates data. The returned list must contain: `dataset`,
#'   `sampling_loops`, `bad_samples`, `testing_data`, and `true_parameters`.
#'
#' @return A list containing the fitted model, dataset, and metadata
#' @keywords internal
#' @export
#'
#' @examples
#' \dontrun{
#' # reproduce_result() requires a custom data_generator_fn
#' # See package vignettes for complete usage examples
#' }
reproduce_result <- function(result, data_generator_fn = NULL) {
  # Check if bayesfam is available
  if (!requireNamespace("bayesfam", quietly = TRUE)) {
    stop("Package 'bayesfam' is required for reproduce_result()")
  }

  family <- bayesfam::brms_family_lookup(
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

  if (is.null(data_generator_fn)) {
    stop(
      "data_generator_fn must be provided. ",
      "Pass a function with signature (config, seed) -> list"
    )
  }

  datagen_result <- data_generator_fn(
    config = data_gen_conf,
    seed = data_gen_conf$seed
  )
  dataset <- datagen_result$dataset
  sampling_loops <- datagen_result$sampling_loops
  bad_samples <- datagen_result$bad_samples
  testing_data <- datagen_result$testing_data

  fit <- stats::update(
    prefit,
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

  list(
    fit = fit,
    dataset = dataset,
    testing_data = datagen_result$testing_data,
    sampling_loops = sampling_loops,
    bad_samples = bad_samples,
    data_gen_conf = data_gen_conf
  )
}
