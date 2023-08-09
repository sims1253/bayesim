# single_model_sim <- function(
#     dataset,
#     fit_configuration,
#     sample_configuration,
#     prefit,
#     seed,
#     debug,
#     result_path,
#     ncores,
#     ...) {
#   fit <- stats::update(prefit,
#     newdata = dataset,
#     formula. = brms::brmsformula(fit_conf$formula),
#     refresh = 0,
#     silent = 2,
#     warmup = stan_pars$warmup,
#     iter = stan_pars$iter,
#     chains = stan_pars$chains,
#     backend = stan_pars$backend,
#     seed = seed,
#     init = stan_pars$init
#   )
#   if (debug == TRUE) {
#     saveRDS(fit, paste0(paste(result_path, "fit", sep = "/"), ".RDS"))
#   }
#
#   all_metric_results <- do.call(
#     metric_list_handler,
#     c(
#       list(
#         fit = fit,
#         data_gen_output = data_gen_output,
#         fit_conf = fit_conf
#       ),
#       list(...)
#     )
#   )
#   if (debug == TRUE) {
#     saveRDS(all_metric_results, paste0(paste(result_path, "metric_results", sep = "/"), ".RDS"))
#   }
#   return(all_metric_results)
# }
#
#
# if (ncores > 1) {
#   if (debug) {
#     cluster <- parallel::makeCluster(ncores,
#       type = cluster_type,
#       outfile = paste(result_path, "cluster_log")
#     )
#   } else {
#     cluster <- parallel::makeCluster(ncores,
#       type = cluster_type
#     )
#   }
#   # Multiprocessing setup
#   doParallel::registerDoParallel(cluster)
#   parallel::clusterEvalQ(cl = cluster, {
#     options(mc.cores = 1)
#   })
#   `%dopar%` <- foreach::`%dopar%`
#
#   # Multiprocessing run
#   final_result <- foreach::foreach(
#     par_seed = seed_list,
#     .packages = "bayesim",
#     .combine = dplyr::bind_rows, .multicombine = TRUE,
#     .maxcombine = length(seed_list),
#     .verbose = debug,
#     .inorder = FALSE
#   ) %dopar% {
#     dataset_sim(
#       data_gen_conf = data_gen_conf,
#       data_gen_fun = data_gen_fun,
#       fit_confs = fit_confs,
#       prefits = prefits,
#       stan_pars = stan_pars,
#       seed = par_seed,
#       debug = debug,
#       result_path = result_path,
#       ...
#     )
#   }
#
#   # Multiprocessing teardown
#   parallel::stopCluster(cluster)
# } else {
#   results <- vector(mode = "list", length = length(seed_list))
#   for (i in seq_along(seed_list)) {
#     results[[i]] <- dataset_sim(
#       data_gen_conf = data_gen_conf,
#       data_gen_fun = data_gen_fun,
#       fit_confs = fit_confs,
#       prefits = prefits,
#       stan_pars = stan_pars,
#       seed = seed_list[[i]],
#       debug = debug,
#       result_path = result_path,
#       ...
#     )
#   }
#   final_result <- dplyr::bind_rows(results)
# }
