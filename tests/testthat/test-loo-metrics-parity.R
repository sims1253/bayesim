# F3 acceptance tests: rmse_loo / r2_loo must match brms::loo_R2 / brms::loo_predict
# on identical draws. Uses a real (tiny) cmdstanr brmsfit.

skip_unless_bayesim_backend()
skip_if_not(requireNamespace("cmdstanr", quietly = TRUE))
skip_if_not(requireNamespace("brms", quietly = TRUE))
skip_if_not(requireNamespace("posterior", quietly = TRUE))
skip_if_not(requireNamespace("loo", quietly = TRUE))

# One shared fixture fit for the whole file (chains = 1, iter = 100).
fit <- suppressWarnings(brms::brm(
  y ~ x,
  data = withr::with_seed(77L, data.frame(y = rnorm(30), x = rnorm(30))),
  family = gaussian(),
  backend = "cmdstanr",
  chains = 1L,
  iter = 100L,
  warmup = 50L,
  seed = 77L,
  silent = 2L,
  refresh = 0L
))

test_that("brms defaults use the fitted rows and explicit newdata stays explicit", {
  train <- fit$data
  train$x[5] <- NA_real_
  dropped <- suppressWarnings(stats::update(
    fit,
    newdata = train,
    recompile = FALSE,
    seed = 77L,
    silent = 2L,
    refresh = 0L
  ))
  result <- list(
    success = TRUE,
    fit = dropped,
    data_bundle = list(train = train)
  )
  fitter <- BrmsFitter()
  expect_equal(log_lik_matrix(fitter, result), brms::log_lik(dropped))
  expect_equal(predict_epred(fitter, result), brms::posterior_epred(dropped))
  withr::local_seed(99L)
  rng <- .Random.seed
  expect_equal(
    predict_fit(fitter, result, seed = 11L)$predicted_samples,
    withr::with_seed(11L, brms::posterior_predict(dropped))
  )
  expect_identical(.Random.seed, rng)
  newdata <- train[1:2, ]
  expect_equal(
    log_lik_matrix(fitter, result, newdata),
    brms::log_lik(dropped, newdata = newdata)
  )
  expect_equal(
    predict_epred(fitter, result, newdata),
    brms::posterior_epred(dropped, newdata = newdata)
  )
})

test_that("brms defaults preserve fitted measurement-error latents", {
  withr::local_package("brms")
  train <- fit$data
  train$x_se <- rep(0.2, nrow(train))
  latent <- suppressWarnings(brms::brm(
    y ~ me(x, x_se),
    data = train,
    backend = "cmdstanr",
    chains = 1L,
    iter = 100L,
    warmup = 50L,
    seed = 77L,
    save_pars = brms::save_pars(latent = TRUE),
    silent = 2L,
    refresh = 0L
  ))
  result <- list(
    success = TRUE,
    fit = latent,
    data_bundle = list(train = train)
  )
  expect_equal(log_lik_matrix(BrmsFitter(), result), brms::log_lik(latent))
  expect_equal(
    predict_epred(BrmsFitter(), result),
    brms::posterior_epred(latent)
  )
})

test_that("brms fit results reject silently dropped training rows", {
  train <- fit$data
  train$x[5] <- NA_real_
  result <- fit_model(
    BrmsFitter(chains = 1L, iter = 100L, warmup = 50L, precompile = FALSE),
    list(train = train, response = "y"),
    list(formula = y ~ x, family = gaussian()),
    seed = 77L,
    task_ctx = list(task_id = "missing-row")
  )
  expect_false(result$success)
  expect_s3_class(result$error, "bayesim_data_error")
  expect_match(conditionMessage(result$error), "training rows")
})

# Build the bayesim context by hand from the fit, mirroring build_loo_context().
ll <- brms::log_lik(fit) # S x N
epred <- brms::posterior_epred(fit) # S x N
y <- as.numeric(brms:::get_y(fit)) # observed response (sorted within brms)
chain_id <- posterior::as_draws_df(fit)$.chain
r_eff <- loo::relative_eff(exp(ll), chain_id = chain_id)
psis_obj <- suppressWarnings(loo::psis(-ll, r_eff = r_eff))
loo_obj <- suppressWarnings(loo::loo(ll))

ctx <- list(
  loo = list(
    elpd = loo_obj$estimates["elpd_loo", "Estimate"],
    p_loo = loo_obj$estimates["p_loo", "Estimate"],
    elpd_se = loo_obj$estimates["elpd_loo", "SE"],
    pareto_k = loo::pareto_k_values(loo_obj)
  ),
  loo_psis = psis_obj,
  loo_psis_ll = ll,
  loo_epred = epred
)
fit_result <- list(fit = fit, success = TRUE)
# Build a data_bundle whose train$response matches the fit's y ordering.
df <- fit$data
df[[1]] <- y # ensure the response column equals brms' get_y() output
resp_name <- all.vars(fit$formula$formula)[1]
data_bundle <- list(
  train = df,
  test = NULL,
  response = resp_name
)

describe("r2_loo_metric parity with brms::loo_R2", {
  it("matches brms::loo_R2 mean within Monte-Carlo tolerance", {
    # brms loo_R2 uses an Exp(1)-weighted draw-wise variance (the posterior R2
    # distribution per Gelman et al. 2018); both brms and our metric draw fresh
    # Exp(1) weights, so exact equality is impossible. Tolerance 0.05 reflects
    # this Monte-Carlo variance on a tiny (100-iteration, 1-chain) fit.
    brms_r2 <- as.numeric(suppressWarnings(brms::loo_R2(fit, summary = TRUE))[
      "R2",
      "Estimate"
    ])
    out <- compute_metric(
      r2_loo_metric(),
      fit_result,
      data_bundle,
      ctx,
      list(task_id = "t1")
    )
    expect_false(is.na(out$value))
    expect_lt(abs(out$value - brms_r2), 0.05)
  })
})

describe("rmse_loo_metric parity with brms loo_predict", {
  it("matches sqrt(mean((y - E_loo_mean)^2))", {
    # brms loo_predict(type='mean') uses posterior_predict; our metric uses epred
    # by default (F3 uses epred for the mean-type LOO prediction). Both target
    # the LOO posterior-predictive mean; for gaussian they coincide in
    # expectation. Compute the E_loo(epred) RMSE directly and compare.
    yloo <- loo::E_loo(epred, psis_obj, log_ratios = -ll, type = "mean")$value
    expected <- sqrt(mean((y - yloo)^2))
    out <- compute_metric(
      rmse_loo_metric(),
      fit_result,
      data_bundle,
      ctx,
      list(task_id = "t1")
    )
    expect_false(is.na(out$value))
    expect_lt(abs(out$value - expected), 1e-6)
    expect_false(is.na(out$pareto_k_max))
  })
})

describe("rmse_loo / r2_loo degradation", {
  it("both NA when PSIS/epred context absent", {
    out_r <- compute_metric(
      rmse_loo_metric(),
      fit_result,
      data_bundle,
      list(loo = ctx$loo),
      list(task_id = "t1")
    )
    out_r2 <- compute_metric(
      r2_loo_metric(),
      fit_result,
      data_bundle,
      list(loo = ctx$loo),
      list(task_id = "t1")
    )
    expect_true(is.na(out_r$value))
    expect_true(is.na(out_r2$value))
  })
})

# loo_fit(BrmsFitter) must use chain-aware r_eff, matching brms::loo().
describe("loo_fit(BrmsFitter) parity with brms::loo", {
  it("matches brms::loo elpd_loo within 1e-6 (chain-aware r_eff)", {
    fitter <- BrmsFitter(chains = 1L, iter = 100L, warmup = 50L, cores = 1L)
    out <- loo_fit(fitter, fit_result)
    brms_loo <- suppressWarnings(brms::loo(fit))
    expect_false(is.na(out$elpd))
    expect_lt(abs(out$elpd - brms_loo$estimates["elpd_loo", "Estimate"]), 1e-6)
    expect_lt(abs(out$p_loo - brms_loo$estimates["p_loo", "Estimate"]), 1e-6)
    expect_lt(abs(out$elpd_se - brms_loo$estimates["elpd_loo", "SE"]), 1e-6)
  })
})

# #73: loo_fit() accepts the train-set log-lik matrix its caller computed and
# returns the r_eff it derived, so build_loo_context() computes each once.
describe("loo_fit(BrmsFitter) with a supplied log_lik matrix (#73)", {
  it("returns the identical summary and the r_eff it used", {
    fitter <- BrmsFitter(chains = 1L, iter = 100L, warmup = 50L, cores = 1L)
    out_passed <- loo_fit(fitter, fit_result, log_lik = ll)
    out_own <- loo_fit(fitter, fit_result)
    expect_equal(out_passed$elpd, out_own$elpd)
    expect_equal(out_passed$p_loo, out_own$p_loo)
    expect_equal(out_passed$elpd_se, out_own$elpd_se)
    expect_equal(out_passed$pareto_k, out_own$pareto_k)
    # The returned r_eff is the chain-aware correction the summary used,
    # matching the fixture's manual relative_eff() computation.
    expect_false(is.null(out_passed$r_eff))
    expect_equal(out_passed$r_eff, r_eff)
  })
})

describe("build_loo_context shares one log-lik and r_eff across summary and PSIS (#73)", {
  it("delivers the same PSIS object as the manual construction", {
    fitter <- BrmsFitter(chains = 1L, iter = 100L, warmup = 50L, cores = 1L)
    ctx <- build_loo_context(fitter, fit_result, need_psis = TRUE)
    expect_equal(ctx$log_lik, ll)
    expect_equal(ctx$loo$elpd, loo_fit(fitter, fit_result)$elpd)
    # The same chain-aware r_eff entered the PSIS smoothing as in the manual
    # fixture (r_eff shifts the tail smoothing, hence pareto_k).
    expect_equal(ctx$psis$diagnostics$pareto_k, psis_obj$diagnostics$pareto_k)
  })
})

describe("build_loo_context PSIS parity", {
  it("matches the PSIS object retained by loo_fit", {
    fitter <- BrmsFitter(chains = 1L, iter = 100L, warmup = 50L, cores = 1L)
    retained <- loo_fit(
      fitter,
      fit_result,
      log_lik = ll,
      save_psis = TRUE
    )$psis_object
    expect_false(is.null(retained))
    ctx <- build_loo_context(fitter, fit_result, need_psis = TRUE)
    # This checks values; test-engine.R separately verifies reuse with a marker.
    expect_identical(ctx$psis, retained)
  })
})
