# Schema-conformance tests for the built-in metric library.
#
# Covers three invariants:
#   1. Degraded (NA) compute() outputs contain only schema-declared fields, so
#      downstream flattening produces stable columns across tasks.
#   2. Field metadata is truthful (pos_prob by_param is a continuous estimate,
#      not a binary outcome; sampler_diagnostics declares all five fields).
#   3. The supports_epred capability flag and validate_fitter() conformance.
#
# All fixtures are synthetic (no brms/cmdstan) so the file stays fast-tier.

.describe_degraded_cases <- function() {
  no_draws_fit <- list(draws = NULL, diagnostics = NULL)
  truthless_bundle <- list(
    train = data.frame(y = rnorm(5), x = rnorm(5)),
    test = NULL,
    response = "y",
    true_params = NULL,
    vars_of_interest = NULL
  )
  list(
    list(pred_rmse_metric(), no_draws_fit, truthless_bundle, list()),
    list(pred_mae_metric(), no_draws_fit, truthless_bundle, list()),
    list(pred_mse_metric(), no_draws_fit, truthless_bundle, list()),
    list(pred_bias_metric(), no_draws_fit, truthless_bundle, list()),
    list(coverage_metric(), no_draws_fit, truthless_bundle, list()),
    list(pos_prob_metric(), no_draws_fit, truthless_bundle, list()),
    list(posterior_summary_metric(), no_draws_fit, truthless_bundle, list()),
    list(sampler_diagnostics_metric(), no_draws_fit, truthless_bundle, list()),
    list(rank_metric(), no_draws_fit, truthless_bundle, list()),
    list(elpd_loo_metric(), no_draws_fit, truthless_bundle, list()),
    list(rmse_loo_metric(), no_draws_fit, truthless_bundle, list()),
    list(r2_loo_metric(), no_draws_fit, truthless_bundle, list()),
    list(elpd_test_metric(), no_draws_fit, truthless_bundle, list()),
    list(r2_test_metric(), no_draws_fit, truthless_bundle, list())
  )
}

describe("degraded metric outputs conform to their schemas", {
  it("degraded field names are a subset of schema field names for every built-in", {
    for (case in .describe_degraded_cases()) {
      metric <- case[[1]]
      out <- compute_metric(
        metric,
        case[[2]],
        case[[3]],
        case[[4]],
        list(task_id = "schema_conformance")
      )
      label <- paste0(metric@name, ": ", paste(names(out), collapse = ", "))
      expect_true(
        all(names(out) %in% names(metric@schema)),
        info = paste("undeclared field in degraded output of", label)
      )
      expect_silent(validate_metric_output(out, metric@name))
      expect_silent(flatten_metric_output(out, metric@name))
    }
  })

  it("coverage degrades to exactly mean/by_param (no stray 'value')", {
    out <- compute_metric(
      coverage_metric(),
      list(draws = NULL, diagnostics = NULL),
      list(response = "y", true_params = NULL),
      list(),
      list(task_id = "t")
    )
    expect_setequal(names(out), c("mean", "by_param"))
    expect_true(is.na(out$mean))
    expect_true(is.na(out$by_param))
  })

  it("pos_prob degrades to exactly mean/by_param", {
    out <- compute_metric(
      pos_prob_metric(),
      list(draws = NULL),
      list(response = "y"),
      list(),
      list(task_id = "t")
    )
    expect_setequal(names(out), c("mean", "by_param"))
    expect_true(is.na(out$mean))
    expect_true(is.na(out$by_param))
  })

  it("posterior_summary degraded output includes the declared median field", {
    out <- compute_metric(
      posterior_summary_metric(),
      list(draws = NULL),
      list(response = "y"),
      list(),
      list(task_id = "t")
    )
    expect_setequal(
      names(out),
      c("mean", "median", "sd", "q_lower", "q_upper")
    )
    expect_true(is.na(out$median))
  })

  it("rank degraded output includes by_param and n_ranks", {
    out <- compute_metric(
      rank_metric(),
      list(draws = NULL),
      list(response = "y", true_params = NULL),
      list(),
      list(task_id = "t")
    )
    expect_setequal(names(out), c("n_draws", "stride", "by_param", "n_ranks"))
    expect_true(is.na(out$by_param))
    expect_true(is.na(out$n_ranks))
  })

  it("degraded outputs pass validate_metric representative execution", {
    # Representative execution checks that every schema-declared field is
    # produced; the NA-degradation paths must satisfy it, not violate it.
    fit_result <- new_fit_result(
      success = TRUE,
      diagnostics = list(),
      timing = list(total = 0)
    )
    bundle <- list(
      train = data.frame(y = rnorm(5), x = rnorm(5)),
      test = NULL,
      response = "y",
      true_params = NULL,
      vars_of_interest = NULL
    )
    expect_silent(validate_metric(
      coverage_metric(),
      fit_result = fit_result,
      data_bundle = bundle
    ))
    expect_silent(validate_metric(
      pos_prob_metric(),
      fit_result = fit_result,
      data_bundle = bundle
    ))
    expect_silent(validate_metric(
      posterior_summary_metric(),
      fit_result = fit_result,
      data_bundle = bundle
    ))
  })

  it("summarize_simulation aggregates all-NA columns to NA without error", {
    df <- data.frame(
      task_id = c("t1", "t2", "t3"),
      status = "success",
      coverage__mean = NA_real_,
      coverage__by_param = NA_real_
    )
    smry <- expect_silent(summarize_simulation(df))
    expect_true(all(is.na(smry$coverage__mean_mean)))
    expect_equal(smry$coverage__mean_n_used, 0L)
  })
})

describe("schema metadata truthfulness", {
  it("pos_prob by_param is an estimate (continuous proportion), not binary", {
    meta <- pos_prob_metric()@schema$by_param
    expect_equal(
      meta,
      list(role = "estimate", aggregation = "mean", mcse = "sd")
    )
  })

  it("coverage by_param remains a binary proportion with nominal level", {
    m <- coverage_metric(prob = 0.9)
    expect_equal(
      m@schema$by_param[c("role", "aggregation", "mcse")],
      list(role = "binary", aggregation = "proportion", mcse = "binomial")
    )
    expect_equal(m@schema$by_param$nominal, 0.9)
  })

  it("sampler_diagnostics declares all five diagnostic fields", {
    m <- sampler_diagnostics_metric()
    expect_setequal(
      names(m@schema),
      c(
        "rhat_max",
        "ess_bulk_min",
        "ess_tail_min",
        "divergent",
        "max_treedepth"
      )
    )
    for (field in c(
      "rhat_max",
      "ess_bulk_min",
      "ess_tail_min",
      "max_treedepth"
    )) {
      expect_equal(
        m@schema[[field]],
        list(role = "diagnostic", aggregation = "mean", mcse = "sd"),
        info = field
      )
    }
    expect_equal(
      m@schema$divergent,
      list(role = "count", aggregation = "none", mcse = "none")
    )
  })
})

describe("supports_epred capability", {
  it("defaults FALSE on the Fitter base class (safe by default)", {
    MinimalFitter <- S7::new_class(
      "MinimalEpredFitter",
      parent = Fitter,
      properties = list(
        name = S7::new_property(S7::class_character, default = "minimal")
      )
    )
    S7::method(fit_model, MinimalFitter) <- function(
      fitter,
      data_bundle,
      fit_spec,
      seed,
      task_ctx
    ) {
      draws <- matrix(rnorm(20), ncol = 2, dimnames = list(NULL, c("a", "b")))
      new_fit_result(
        success = TRUE,
        draws = draws,
        diagnostics = list(),
        timing = list(total = 0.1)
      )
    }
    S7::method(extract_draws, MinimalFitter) <- function(
      fitter,
      fit_result,
      variables = NULL
    ) {
      fit_result$draws
    }
    fitter <- MinimalFitter()
    expect_false(fitter@supports_epred)
    expect_null(predict_epred(fitter, list()))
    expect_silent(validate_fitter(fitter))
  })

  it("LinearRegressionFitter declares and provides epred (S x N)", {
    fitter <- LinearRegressionFitter(n_draws = 40L)
    expect_true(fitter@supports_epred)

    bundle <- list(
      train = data.frame(y = rnorm(12), x = rnorm(12)),
      test = NULL,
      response = "y"
    )
    fit_result <- fit_model(
      fitter,
      bundle,
      list(model = "lm"),
      seed = 11L,
      task_ctx = list(task_id = "epred_test")
    )
    epred <- predict_epred(fitter, fit_result)
    expect_true(is.matrix(epred))
    expect_true(is.numeric(epred))
    expect_equal(dim(epred), c(40L, 12L)) # S x N: draws x observations
  })

  it("validate_fitter smoke test exercises epred when declared TRUE", {
    expect_silent(validate_fitter(
      LinearRegressionFitter(n_draws = 25L),
      smoke_test = TRUE
    ))
  })

  it("validate_fitter rejects supports_epred = TRUE without an override", {
    BadEpredFitter <- S7::new_class(
      "BadEpredFitter",
      parent = Fitter,
      properties = list(
        name = S7::new_property(S7::class_character, default = "bad_epred"),
        supports_epred = S7::new_property(S7::class_logical, default = TRUE)
      )
    )
    S7::method(fit_model, BadEpredFitter) <- function(
      fitter,
      data_bundle,
      fit_spec,
      seed,
      task_ctx
    ) {
      draws <- matrix(rnorm(20), ncol = 2, dimnames = list(NULL, c("a", "b")))
      new_fit_result(
        success = TRUE,
        draws = draws,
        diagnostics = list(),
        timing = list(total = 0.1)
      )
    }
    S7::method(extract_draws, BadEpredFitter) <- function(
      fitter,
      fit_result,
      variables = NULL
    ) {
      fit_result$draws
    }
    expect_error(
      validate_fitter(BadEpredFitter()),
      "supports_epred"
    )
  })

  it("CmdStanFitter derives supports_epred from the epred GQ argument", {
    stan_data_fn <- function(data_bundle, fit_spec) list()
    without <- CmdStanFitter(
      stan_code = "parameters { real y; } model { y ~ normal(0, 1); }",
      stan_data = stan_data_fn
    )
    expect_false(without@supports_epred)
    with_epred <- CmdStanFitter(
      stan_code = "parameters { real y; } model { y ~ normal(0, 1); }",
      stan_data = stan_data_fn,
      epred = "mu"
    )
    expect_true(with_epred@supports_epred)
  })

  it("CmdStanFitter predict_fit returns NULL (seam warns), not an error", {
    fitter <- CmdStanFitter(
      stan_code = "parameters { real y; } model { y ~ normal(0, 1); }",
      stan_data = function(data_bundle, fit_spec) list()
    )
    expect_null(predict_fit(fitter, list()))
  })
})

describe("error condition classes", {
  it("bayesim_contract_error is exported with the documented class chain", {
    expect_true("bayesim_contract_error" %in% getNamespaceExports("bayesim"))
    err <- bayesim_contract_error("Contract violation")
    expect_s3_class(err, "bayesim_contract_error")
    expect_s3_class(err, "bayesim_error")
    expect_s3_class(err, "error")
    expect_true(is_fatal_error(err))
  })

  it("bayesim_validation_error is a bayesim_contract_error subclass", {
    err <- bayesim_validation_error("Validation failed")
    expect_s3_class(err, "bayesim_validation_error")
    expect_s3_class(err, "bayesim_contract_error")
    expect_s3_class(err, "bayesim_error")
  })

  it("validate_metric failures are catchable as bayesim_contract_error", {
    draws <- matrix(rnorm(20), ncol = 1, dimnames = list(NULL, "b_x"))
    fit_result <- new_fit_result(
      success = TRUE,
      draws = draws,
      diagnostics = list(),
      timing = list(total = 0.1)
    )
    bundle <- list(
      train = data.frame(y = rnorm(5), x = rnorm(5)),
      test = NULL,
      response = "y",
      true_params = c(nonexistent = 1),
      vars_of_interest = "nonexistent"
    )
    expect_error(
      validate_metric(
        coverage_metric(),
        fit_result = fit_result,
        data_bundle = bundle
      ),
      class = "bayesim_contract_error"
    )
  })
})
