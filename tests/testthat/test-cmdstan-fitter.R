# D2: CmdStanFitter — user-supplied Stan programs. cmdstan-gated.
skip_on_cran()
skip_if_not(requireNamespace("cmdstanr", quietly = TRUE))
skip_if_not(nzchar(Sys.which("cmdstan")) || !is.null(cmdstanr::cmdstan_version()))

stan_file <- system.file("stan", "linear_regression.stan", package = "bayesim")

# stan_data: build the design matrix from a data bundle (response + predictors).
stan_data_fn <- function(data_bundle, fit_spec) {
  data <- data_bundle$train
  response <- data_bundle$response
  predictors <- setdiff(names(data), response)
  X <- cbind(1, as.matrix(data[, predictors, drop = FALSE]))
  list(
    N = nrow(data),
    K = ncol(X) - 1L,
    X = X,
    y = data[[response]]
  )
}

.gen <- function(data_spec, task_ctx) {
  n <- data_spec$n %||% 40L
  x <- stats::rnorm(n)
  y <- 1 + 0.5 * x + stats::rnorm(n)
  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(beta_intercept = 1, beta_x = 0.5, sigma = 1),
    vars_of_interest = c("beta_intercept", "beta_x", "sigma"),
    references = c(beta_intercept = 0, beta_x = 0, sigma = 1),
    meta = list()
  )
}

describe("CmdStanFitter fit + draws + diagnostics", {
  it("fits and recovers parameters", {
    skip_if_not(file.exists(stan_file), "stan file not shipped")
    fitter <- CmdStanFitter(
      stan_file = stan_file,
      stan_data = stan_data_fn,
      log_lik = "log_lik",
      epred = "mu",
      chains = 1L, iter_warmup = 100L, iter_sampling = 100L
    )
    data_bundle <- .gen(list(n = 40L), list(task_id = "t"))
    res <- fit_model(fitter, data_bundle, list(), seed = 1L,
                     task_ctx = list(task_id = "t"))
    expect_true(res$success)

    draws <- extract_draws(fitter, res)
    expect_true(is.matrix(draws))
    expect_true("sigma" %in% colnames(draws))

    # Diagnostics present.
    diag <- fit_diagnostics(fitter, res)
    expect_true(all(c("rhat_max", "ess_bulk_min", "divergent") %in% names(diag)))
    expect_false(is.na(diag$rhat_max))

    # log_lik is S x N.
    ll <- log_lik_matrix(fitter, res)
    expect_equal(dim(ll), c(100L, 40L))

    # epred is S x N.
    epred <- predict_epred(fitter, res)
    expect_equal(dim(epred), c(100L, 40L))

    # loo works.
    loo <- loo_fit(fitter, res)
    expect_false(is.na(loo$elpd))
  })

  it("errors when log_lik is needed but not declared", {
    skip_if_not(file.exists(stan_file), "stan file not shipped")
    fitter <- CmdStanFitter(
      stan_file = stan_file,
      stan_data = stan_data_fn,
      chains = 1L, iter_warmup = 50L, iter_sampling = 50L
    )
    data_bundle <- .gen(list(n = 30L), list(task_id = "t"))
    res <- fit_model(fitter, data_bundle, list(), seed = 1L,
                     task_ctx = list(task_id = "t"))
    expect_error(log_lik_matrix(fitter, res), class = "bayesim_config_error")
  })

  it("runs under workers = 2 (matches sequential summary)", {
    skip_if_not(file.exists(stan_file), "stan file not shipped")
    fitter <- CmdStanFitter(
      stan_file = stan_file,
      stan_data = stan_data_fn,
      log_lik = "log_lik",
      chains = 1L, iter_warmup = 50L, iter_sampling = 50L
    )
    config <- simulation_config(
      data_grid = data.frame(n = 30L),
      fit_grid = data.frame(model = "stan"),
      data_generator = .gen,
      fitter = fitter,
      metrics = list(),
      n_replicates = 2L,
      seed = 7L
    )
    seq_res <- run_simulation(config, resume = "never", progress = FALSE)
    par_res <- run_simulation(config, resume = "never", progress = FALSE, workers = 2L)
    expect_false(mirai::daemons_set())
    expect_equal(nrow(seq_res$summary), nrow(par_res$summary))
    expect_true(all(par_res$summary$status == "success"))
  })
})
