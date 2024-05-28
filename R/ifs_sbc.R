#' Full inverse forward sampling supported SBC
#'
#' @param fit brmsfit object that is used to generate SBC datasets and fit SBC models
#' @param n_sims Number of SBC simulations/datasets.
#' @param ppred_data_gen Function that generates the data for the newdata argument of posterior_predict to generate SBC datasets.
#' @param precon_sample The dataset that was used to precondition the fit.
#' @param lb Lower bound for the outcome. Is used to resample.
#' @param ub Upper bound for the outcome. Is used to resample.
#' @param truncate If TRUE, will truncate outcome values, if FALSE, will resample instead.
#' @param post_warmup_samples Is passed to the sbc model updates and used to
#' make sure resampling only uses valid samples. Increase in case of numerical instability.
#' @param ... further parameters. Currently only passed to ppred_data_gen.
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
                    truncate = FALSE,
                    post_warmup_samples = 1000,
                    ...) {
  model_parameters <- c(posterior::variables(fit), "loglik")
  model_parameters <- model_parameters[!model_parameters %in% c(
    "lp__",
    "lprior"
  )]

  # Prepare rank summary DF
  ranks <- as.data.frame(
    matrix(nrow = n_sims, ncol = length(model_parameters))
  )
  colnames(ranks) <- model_parameters

  # prepare draw indices
  index_list <- sample(1:ndraws(fit), size = ndraws(fit), replace = FALSE)

  results <- future.apply::future_lapply(
    index_list[1:n_sims],
    sbc_sim,
    fit = fit,
    ppred_data_gen = ppred_data_gen,
    precon_sample = precon_sample,
    lb = lb,
    ub = ub,
    truncate = truncate,
    post_warmup_samples = post_warmup_samples,
    ...,
    future.seed = TRUE,
    future.chunk.size = floor(SBC::default_chunk_size(n_sims)),
    future.packages = c("bayesim")
  )

  if (length(
    bad_indices <- which(
      sapply(results, function(x) any(is.na(x)))
    )
  ) > (n_sims / (post_warmup_samples / n_sims))) {
    warning("Rate of numerical instability is too high and likely won't be able
         to be resolved via resampling.")
  }

  # If we want to resample, we do so until we are out of the posterior samples.
  if (!truncate) {
    i <- n_sims
    while (length(
      bad_indices <- which(
        sapply(results, function(x) any(is.na(x)))
      )
    ) > 0
    ) {
      if (i + length(bad_indices) > post_warmup_samples) {
        warning("Resampling ran out of samples which indicates a high degree of
           numerical instability in the data generating process. Retry with
           more samples, increase the size of the preconditioning dataset or
           use narrower priors.")
        results[bad_indices[1:(post_warmup_samples - i)]] <- future.apply::future_lapply(
          index_list[(i + 1):post_warmup_samples],
          sbc_sim,
          fit = fit,
          ppred_data_gen = ppred_data_gen,
          precon_sample = precon_sample,
          lb = lb,
          ub = ub,
          truncate = truncate,
          post_warmup_samples = post_warmup_samples,
          ...,
          future.seed = TRUE,
          future.chunk.size = floor(SBC::default_chunk_size(n_sims)),
          future.packages = c("bayesim")
        )
        results <- results[!which(
          sapply(results, function(x) any(is.na(x)))
        )]
        break # there is proably a cleaner way to do this...
      } else {
        results[bad_indices] <- future.apply::future_lapply(
          index_list[(i + 1):(i + length(bad_indices))],
          sbc_sim,
          fit = fit,
          ppred_data_gen = ppred_data_gen,
          precon_sample = precon_sample,
          lb = lb,
          ub = ub,
          truncate = truncate,
          post_warmup_samples = post_warmup_samples,
          ...,
          future.seed = TRUE,
          future.chunk.size = floor(SBC::default_chunk_size(n_sims)),
          future.packages = c("bayesim")
        )
        i <- i + length(bad_indices)
      }
    }
  }
  # resample until we found enough working cases or ran out of samples


  # Clean up the results to work with SBC plotting functions
  for (i in seq_along(results)) {
    ranks[i, ] <- results[[i]]
  }
  ranks_df <- ranks %>%
    mutate(sim_id = seq_len(n())) %>%
    pivot_longer(c(-sim_id), names_to = "variable", values_to = "rank")
  return(list("ranks_df" = ranks_df, diagnostics = list("resample" = i)))
}

sbc_sim <- function(index, fit, ppred_data_gen, precon_sample, lb, ub, ...,
                    truncate = FALSE, post_warmup_samples = 1000) {
  gen_dataset <- brms_full_ppred(
    fit = fit,
    draws = index,
    newdata = ppred_data_gen(fit, ...)
  )[[index]]

  responses <- brms_response_sequence(fit)
  if (!truncate) {
    # If any of the response variables has a single entry that is outside the
    # boundaries, we return a frame of NAs to force resampling.
    if (length(
      which(
        if (is.infinite(lb)) {
          is.infinite(unlist(gen_dataset[unlist(responses)]))
        } else {{
          unlist(gen_dataset[unlist(responses)] <= lb)
        } |
          if (is.infinite(ub)) {
            is.infinite(unlist(gen_dataset[unlist(responses)]))
          } else {
            unlist(gen_dataset[unlist(responses)] >= ub)
          }}
      )
    ) > 0) {
      model_parameters <- c(posterior::variables(fit), "loglik")
      model_parameters <- model_parameters[!model_parameters %in% c(
        "lp__",
        "lprior"
      )]
      tmp <- numeric(length = length(model_parameters))
      tmp <- sapply(tmp, function(x) NA)
      names(tmp) <- model_parameters
      return(tmp)
    }
  } else {
    # Trunkates the responses to the given boundaries.
    for (resp in responses)
    {
      gen_dataset[[resp]] <- ifelse(gen_dataset[[resp]] < lb, lb, gen_dataset[[resp]])
      gen_dataset[[resp]] <- ifelse(gen_dataset[[resp]] > ub, ub, gen_dataset[[resp]])
    }
  }


  # Collect true parameters for SBC rank comparison
  true_pars <- as.data.frame(posterior::as_draws_matrix(fit)[index, ])
  true_ll <- sum(brms::log_lik(fit,
    gen_dataset,
    draw_ids = index
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
    cores = 1,
    refresh = 0,
    silent = 2,
    iter = 2 * post_warmup_samples,
    sample_prior = "no" # this is important so we don't accidentally only sample the prior
  )
  fit_pars <- as.data.frame(as_draws_matrix(sbc_fit))
  fit_pars <- fit_pars[!names(fit_pars) %in% c("lprior", "lp__")]
  fit_ll <- rowSums(log_lik(sbc_fit, gen_dataset))

  # Calculate the rank statistics by comparing the true parameters with the
  # draws from the sbc fit.
  model_parameters <- c(posterior::variables(fit), "loglik")
  model_parameters <- model_parameters[!model_parameters %in% c(
    "lp__",
    "lprior"
  )]
  tmp <- vector(mode = "numeric", length = length(model_parameters))
  names(tmp) <- model_parameters
  for (name in colnames(fit_pars)) {
    tmp[[name]] <- sum(fit_pars[[name]] < true_pars[[name]])
  }
  tmp[["loglik"]] <- sum(fit_ll < true_ll)
  return(tmp)
}
