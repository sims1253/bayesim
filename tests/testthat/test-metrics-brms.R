# F2 acceptance tests: metrics must produce non-NA output when
# vars_of_interest holds cleaned (non-b_-prefixed) names but the draws
# matrix carries the b_ prefix (the real BrmsFitter condition).

describe("resolve_draw_columns()", {
  it("maps plain names unchanged when present", {
    m <- resolve_draw_columns(c("alpha", "beta"), c("alpha", "beta", "sigma"))
    expect_equal(unname(m), c("alpha", "beta"))
    expect_equal(names(m), c("alpha", "beta"))
  })

  it("maps cleaned names to b_-prefixed columns", {
    m <- resolve_draw_columns(c("x", "Intercept"),
                              c("b_x", "b_Intercept", "sigma"))
    expect_equal(unname(m), c("b_x", "b_Intercept"))
    expect_equal(names(m), c("x", "Intercept"))
  })

  it("errors with bayesim_config_error when a var is absent", {
    expect_error(
      resolve_draw_columns(c("x", "nonexistent"), c("b_x", "sigma")),
      class = "bayesim_config_error"
    )
  })

  it("handles empty input", {
    m <- resolve_draw_columns(character(0), c("b_x"))
    expect_equal(length(m), 0L)
  })
})

describe("metrics against a real (tiny) brmsfit", {
  skip_on_cran()
  skip_if_not(requireNamespace("cmdstanr", quietly = TRUE))
  skip_if_not(requireNamespace("brms", quietly = TRUE))
  skip_if_not(requireNamespace("posterior", quietly = TRUE))

  fit <- brms::brm(
    y ~ x,
    data = data.frame(y = rnorm(30), x = rnorm(30)),
    family = gaussian(),
    backend = "cmdstanr",
    chains = 1L, iter = 100L, warmup = 50L,
    silent = 2L, refresh = 0L
  )

  fit_result <- list(draws = brms::as_draws_matrix(fit), diagnostics = NULL)
  data_bundle <- list(
    train = data.frame(y = rnorm(10), x = rnorm(10)),
    test = NULL,
    response = "y",
    vars_of_interest = c("x", "Intercept", "sigma"),  # cleaned names — the F2 condition
    true_params = c(x = 0.5, Intercept = 0.5, sigma = 1)
  )
  ctx <- list(task_id = "t1")

  it("coverage_metric gives 0/1 (not NA) for x/Intercept/sigma", {
    out <- compute(coverage_metric(), fit_result, data_bundle, list(), ctx)
    for (v in c("x", "Intercept", "sigma")) {
      expect_true(v %in% names(out$by_param),
                  info = paste("missing", v))
      expect_false(is.na(out$by_param[[v]]),
                   info = paste("NA for", v))
      expect_true(out$by_param[[v]] %in% c(0, 1),
                  info = paste("non-binary for", v))
    }
  })

  it("posterior_mean_metric returns finite values for x/Intercept/sigma", {
    out <- compute(posterior_mean_metric(), fit_result, data_bundle, list(), ctx)
    for (v in c("x", "Intercept", "sigma")) {
      expect_true(v %in% names(out), info = paste("missing", v))
      expect_true(is.finite(out[[v]]), info = paste("non-finite for", v))
    }
  })

  it("pos_prob_metric returns in-range non-NA probabilities", {
    out <- compute(pos_prob_metric(), fit_result, data_bundle, list(), ctx)
    for (v in c("x", "Intercept", "sigma")) {
      expect_true(v %in% names(out$by_param), info = paste("missing", v))
      expect_false(is.na(out$by_param[[v]]), info = paste("NA for", v))
      expect_gte(out$by_param[[v]], 0)
      expect_lte(out$by_param[[v]], 1)
    }
  })

  it("posterior_summary_metric returns non-NA mean/sd/q_lower/q_upper", {
    out <- compute(posterior_summary_metric(), fit_result, data_bundle, list(), ctx)
    for (field in c("mean", "sd", "q_lower", "q_upper")) {
      for (v in c("x", "Intercept", "sigma")) {
        expect_true(v %in% names(out[[field]]),
                    info = paste("missing", field, v))
        expect_true(is.finite(out[[field]][[v]]),
                    info = paste("non-finite", field, v))
      }
    }
  })

  it("rank_metric returns integer non-NA ranks", {
    out <- compute(rank_metric(), fit_result, data_bundle, list(), ctx)
    for (v in c("x", "Intercept", "sigma")) {
      expect_true(v %in% names(out$by_param), info = paste("missing", v))
      expect_false(is.na(out$by_param[[v]]), info = paste("NA for", v))
      expect_true(is.integer(out$by_param) || is.numeric(out$by_param))
    }
  })

  it("coverage_metric errors when a var is genuinely absent", {
    bad_bundle <- data_bundle
    bad_bundle$vars_of_interest <- c("x", "nonexistent")
    bad_bundle$true_params <- c(x = 0.5, nonexistent = 0.5)
    expect_error(
      compute(coverage_metric(), fit_result, bad_bundle, list(), ctx),
      class = "bayesim_config_error"
    )
  })
})
