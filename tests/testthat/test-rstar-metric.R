# Tests for rstar_metric (ROADMAP I2).
#
# rstar_metric wraps posterior::rstar(), which needs per-chain posterior draws
# and the caret + (ranger|gbm) classifier stack. We exercise two paths:
#   1. A fitter with no chain info (LinearRegressionFitter): must return NA.
#   2. A real (tiny) brmsfit with >= 2 chains: must return a finite value in
#      ~[1, 1.2] (gated on cmdstanr + caret + ranger availability).

describe("rstar_metric construction", {
  it("is an RstarMetric with default name and scalar summary", {
    m <- rstar_metric()
    expect_s7_class(m, RstarMetric)
    expect_equal(m@name, "rstar")
    expect_equal(m@needs, character())
    expect_false(m@required)
    expect_false(m@uncertainty)
  })

  it("accepts uncertainty and method overrides", {
    m <- rstar_metric(uncertainty = TRUE, method = "gbm")
    expect_true(m@uncertainty)
    expect_equal(m@method, "gbm")
  })
})

describe("rstar_metric against a chain-less fitter (LinearRegressionFitter)", {
  skip_if_not(requireNamespace("posterior", quietly = TRUE))

  # LinearRegressionFitter produces a fit_result whose $fit is a plain list
  # with no chain info: rstar_metric must degrade to NA with a warning.
  it("returns NA_real_ when no per-chain draws are available", {
    fit_result <- list(
      fit = list(fitter = "linear_regression"), # no draws() method, no $fit
      draws = matrix(
        rnorm(200),
        ncol = 2,
        dimnames = list(NULL, c("b_x", "sigma"))
      )
    )
    data_bundle <- list(response = "y")
    expect_warning(
      out <- compute_metric(
        rstar_metric(),
        fit_result,
        data_bundle,
        list(),
        list(task_id = "t")
      ),
      regexp = "no chain info|>= 2 chains"
    )
    expect_true(is.na(out$value))
  })
})

describe("rstar_metric against a real brmsfit", {
  skip_unless_bayesim_backend()
  skip_if_not(requireNamespace("cmdstanr", quietly = TRUE))
  skip_if_not(requireNamespace("brms", quietly = TRUE))
  skip_if_not(requireNamespace("posterior", quietly = TRUE))
  skip_if_not(requireNamespace("caret", quietly = TRUE))
  skip_if_not(requireNamespace("ranger", quietly = TRUE))
  skip_if_not(requireNamespace("randomForest", quietly = TRUE))
  skip_if_not(
    nzchar(Sys.which("cmdstan")) ||
      !is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL))
  )

  fit <- brms::brm(
    y ~ x,
    data = data.frame(y = rnorm(40), x = rnorm(40)),
    family = gaussian(),
    backend = "cmdstanr",
    chains = 2L,
    iter = 200L,
    warmup = 100L,
    silent = 2L,
    refresh = 0L
  )

  fit_result <- list(
    fit = fit,
    draws = brms::as_draws_matrix(fit)
  )
  data_bundle <- list(response = "y")
  ctx <- list(task_id = "t1")

  it("returns a finite numeric R* value in a reasonable range", {
    out <- compute_metric(
      rstar_metric(),
      fit_result,
      data_bundle,
      ctx,
      list(task_id = "t1")
    )
    expect_true(
      is.numeric(out$value),
      info = paste("value:", paste(head(out$value), collapse = ","))
    )
    expect_true(
      is.finite(out$value),
      info = paste("non-finite value:", out$value)
    )
    # R* is 2x classifier accuracy; ~1 for converged chains, but sampling
    # noise in the classifier can push it slightly below 1.
    expect_gte(out$value, 0.8)
    expect_lt(out$value, 1.5)
  })

  it("averages the uncertainty vector when uncertainty = TRUE", {
    set.seed(123)
    out <- compute_metric(
      rstar_metric(uncertainty = TRUE),
      fit_result,
      data_bundle,
      ctx,
      list(task_id = "t1")
    )
    expect_true(is.numeric(out$value))
    expect_true(is.finite(out$value))
    expect_gte(out$value, 0.8)
    expect_lt(out$value, 1.5)
  })
})
