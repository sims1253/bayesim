#' Full inverse forward sampling supported SBC
#'
#' @param fit
#' @param n_sims
#' @param ppred_data_gen
#' @param precon_sample
#' @param lb
#' @param ub

#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ifs_SBC <- function(fit,
                    n_sims,
                    ppred_data_gen,
                    precon_sample = NULL,
                    lb = NULL,
                    ub = NULL,
                    ...) {
  model_parameters <- c(posterior::variables(fit), "loglik")

  # Prepare rank summary DF
  ranks <- as.data.frame(
    matrix(nrow = n_sims, ncol = length(model_parameters))
  )
  colnames(ranks) <- model_parameters

  # prepare draw indices
  index_list <- sample(1:ndraws(fit), size = n_sims, replace = FALSE)

  results <- future.apply::future_lapply(
    index_list,
    sbc_sim,
    fit = fit,
    ppred_data_gen = ppred_data_gen,
    precon_sample = precon_sample,
    lb = lb,
    ub = ub,
    ...,
    future.seed = TRUE,
    future.chunk.size = floor(SBC::default_chunk_size(n_sims)),
    future.packages = c("bayesim")
  )

  # Clean up the results to work with SBC plotting functions
  for (i in seq_along(results)) {
    ranks[i, ] <- results[[i]]
  }
  ranks_df <- ranks %>%
    mutate(sim_id = seq_len(n())) %>%
    pivot_longer(c(-sim_id), names_to = "variable", values_to = "rank")

  return(ranks_df)
}


sbc_sim <- function(index, fit, ppred_data_gen, precon_sample, lb, ub, ...) {
  options(mc.cores = 1)

  gen_dataset <- brms_full_ppred(
    fit = fit,
    draws = index,
    newdata = ppred_data_gen(fit, ...)
  )[[index]]

  responses <- brms_response_sequence(fit)
  resample_counter <- 0

  while (length(
    bad_indices <- which(
      if (is.na(lb)) {
        is.infinite(unlist(gen_dataset[unlist(responses)]))
      } else {
        {
          unlist(gen_dataset[unlist(responses)] <= lb)
        } |
          if (is.na(ub)) {
            is.infinite(unlist(gen_dataset[unlist(responses)]))
          } else {
            unlist(gen_dataset[unlist(responses)] >= ub)
          }
      }
    )
  ) > 0
  ) {
    if (resample_counter > 100) stop("Couldn't sample fitting values after 100
                                       resample tries.")
    resample_counter <- resample_counter + 1

    gen_dataset[bad_indices, ] <- brms_full_ppred(
      fit = fit,
      draws = index,
      newdata = ppred_data_gen(fit, ...)[bad_indices, ]
    )[[index]]
  }

  # Collect true parameters for SBC rank comparison
  true_pars <- as.data.frame(posterior::as_draws_matrix(fit)[index, ])
  true_ll <- sum(brms::log_lik(fit,
    gen_dataset,
    draw_ids = index,
    allow_new_levels = TRUE
  ))

  # Preserve the generated dataset for the loglik later.
  full_data <- gen_dataset
  # Add the precon sample to the dataset if provided
  if (!is.null(precon_sample)) {
    full_data <- rbind(precon_sample, full_data)
  }

  # Fit model to sbc dataset and extract the poserior draws and log likelihood
  sbc_fit <- update(fit,
    newdata = full_data,
    chains = 1,
    init_r = 0.1,
    refresh = 0,
    silent = 2
  )
  fit_pars <- as.data.frame(as_draws_matrix(sbc_fit))
  fit_ll <- rowSums(log_lik(sbc_fit, gen_dataset))

  # Calculate the rank statistics by comparing the true parameters with the
  # draws from the sbc fit.
  model_parameters <- c(posterior::variables(fit), "loglik")
  tmp <- vector(mode = "numeric", length = length(model_parameters))
  names(tmp) <- model_parameters
  for (name in colnames(fit_pars)) {
    tmp[[name]] <- sum(fit_pars[[name]] < true_pars[[name]])
  }
  tmp[["loglik"]] <- sum(fit_ll < true_ll)
  # tmp[["resample_counter"]] <- resample_counter
  return(tmp)
}
