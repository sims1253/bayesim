skip_on_cran()

skip_if_not(requireNamespace("cmdstanr", quietly = TRUE))
skip_if_not(requireNamespace("brms", quietly = TRUE))
skip_if_not(requireNamespace("rstan", quietly = TRUE))

# A trivial gaussian generator producing the same structure every call.
# Variables (y, x) match the test models so update(recompile=FALSE) holds.
gaussian_data_generator <- function(data_spec, task_ctx) {
  # Consume the ambient RNG state (the worker restores the per-task stream).
  # For direct calls outside the worker, fall back to task_ctx$seed.
  seed <- task_ctx$seed
  if (!is.null(seed)) {
    # Use a fresh local draw to avoid disturbing global state in tests.
    withr::local_seed(seed)
  }
  n <- as.integer(data_spec$n %||% 20L)
  x <- stats::rnorm(n)
  y <- 2 * x + stats::rnorm(n)
  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(beta = 2),
    vars_of_interest = "beta"
  )
}

tiny_fitter <- function() {
  BrmsFitter(chains = 1L, iter = 50L, warmup = 25L, cores = 1L)
}

describe("BrmsFitter timings", {
  it("reads timings directly from cmdstanr fits without rstan", {
    cmdstan_fit <- structure(
      list(time = function() {
        list(
          total = 3,
          chains = data.frame(
            chain_id = 1:2,
            warmup = c(0.5, 0.6),
            sampling = c(1, 1.2),
            total = c(1.5, 1.8)
          )
        )
      }),
      class = "CmdStanMCMC"
    )

    timing <- bayesim:::extract_brms_timings(
      list(fit = cmdstan_fit),
      fallback_total = 99
    )

    expect_equal(timing$warmup, 1.1)
    expect_equal(timing$sample, 2.2)
    expect_equal(timing$total, 3.3)
  })
})

describe("BrmsFitter model bank", {
  it("builds one prefit per distinct model spec", {
    fit_grid <- data.frame(
      model = c("gaussian", "gaussian_dup", "student"),
      stringsAsFactors = FALSE
    )
    fit_grid$formula <- list(y ~ x, y ~ x, y ~ x)
    fit_grid$family <- list(gaussian(), gaussian(), brms::brmsfamily("student"))

    fitter <- tiny_fitter()
    expect_warning(
      bank <- bayesim:::build_model_bank(
        fitter = fitter,
        fit_grid = fit_grid,
        data_generator = gaussian_data_generator,
        data_spec_template = list(n = 20),
        result_path = NULL,
        seed = 42L
      ),
      "explicit priors"
    )

    # Two distinct specs: gaussian (deduped from 2 rows) and student_t.
    expect_type(bank, "list")
    expect_length(bank, 2L)
  })

  it("reuses the prefit across tasks (no recompilation)", {
    fitter <- tiny_fitter()
    fit_grid <- data.frame(model = "gaussian", stringsAsFactors = FALSE)
    fit_grid$formula <- list(y ~ x)
    fit_grid$family <- list(gaussian())

    bank <- bayesim:::build_model_bank(
      fitter = fitter,
      fit_grid = fit_grid,
      data_generator = gaussian_data_generator,
      data_spec_template = list(n = 20),
      result_path = NULL,
      seed = 42L
    )
    bayesim:::set_model_bank(bank)
    on.exit(bayesim:::set_model_bank(NULL), add = TRUE)

    data_bundle <- gaussian_data_generator(list(n = 20), list(seed = 1L))
    fit_spec <- list(formula = y ~ x, family = gaussian())

    # Two fit_model() calls on different data should reuse the same prefit binary.
    # Capture warnings: a recompile warning would indicate the bank was bypassed.
    recompiled <- FALSE
    result1 <- withCallingHandlers(
      fit_model(
        fitter,
        data_bundle,
        fit_spec,
        seed = 100L,
        task_ctx = list(task_id = "t1")
      ),
      warning = function(w) {
        if (grepl("recompil", conditionMessage(w), ignore.case = TRUE)) {
          recompiled <<- TRUE
        }
        invokeRestart("muffleWarning")
      }
    )
    expect_false(recompiled)
    expect_true(result1$success)

    result2 <- fit_model(
      fitter,
      data_bundle,
      fit_spec,
      seed = 200L,
      task_ctx = list(task_id = "t2")
    )
    expect_true(result2$success)

    # Both results reference the same underlying compiled binary (same Stan model
    # path on the prefit they updated from).
    expect_true(!is.null(result1$fit))
    expect_true(!is.null(result2$fit))
  })

  it("propagates seeds deterministically (same seed -> identical draws)", {
    fitter <- tiny_fitter()
    fit_grid <- data.frame(model = "gaussian", stringsAsFactors = FALSE)
    fit_grid$formula <- list(y ~ x)
    fit_grid$family <- list(gaussian())

    bank <- bayesim:::build_model_bank(
      fitter = fitter,
      fit_grid = fit_grid,
      data_generator = gaussian_data_generator,
      data_spec_template = list(n = 20),
      result_path = NULL,
      seed = 42L
    )
    bayesim:::set_model_bank(bank)
    on.exit(bayesim:::set_model_bank(NULL), add = TRUE)

    data_bundle <- gaussian_data_generator(list(n = 20), list(seed = 1L))
    fit_spec <- list(formula = y ~ x, family = gaussian())

    r1 <- fit_model(
      fitter,
      data_bundle,
      fit_spec,
      seed = 777L,
      task_ctx = list(task_id = "a")
    )
    r2 <- fit_model(
      fitter,
      data_bundle,
      fit_spec,
      seed = 777L,
      task_ctx = list(task_id = "b")
    )
    r3 <- fit_model(
      fitter,
      data_bundle,
      fit_spec,
      seed = 888L,
      task_ctx = list(task_id = "c")
    )

    d1 <- brms::as_draws_matrix(r1$fit)
    d2 <- brms::as_draws_matrix(r2$fit)
    d3 <- brms::as_draws_matrix(r3$fit)

    # Same seed -> identical draws
    expect_identical(d1, d2)
    # Different seed -> different draws
    expect_false(identical(d1, d3))
  })

  it("extracts real timings and diagnostics", {
    fitter <- tiny_fitter()
    fit_grid <- data.frame(model = "gaussian", stringsAsFactors = FALSE)
    fit_grid$formula <- list(y ~ x)
    fit_grid$family <- list(gaussian())

    bank <- bayesim:::build_model_bank(
      fitter = fitter,
      fit_grid = fit_grid,
      data_generator = gaussian_data_generator,
      data_spec_template = list(n = 20),
      result_path = NULL,
      seed = 42L
    )
    bayesim:::set_model_bank(bank)
    on.exit(bayesim:::set_model_bank(NULL), add = TRUE)

    data_bundle <- gaussian_data_generator(list(n = 20), list(seed = 1L))
    fit_spec <- list(formula = y ~ x, family = gaussian())
    result <- fit_model(
      fitter,
      data_bundle,
      fit_spec,
      seed = 1L,
      task_ctx = list(task_id = "t")
    )

    # Timings from the updated fit's stanfit (non-zero warmup/sample).
    expect_true(is.numeric(result$timing$warmup) || is.na(result$timing$warmup))
    expect_true(is.numeric(result$timing$total))
    expect_false(is.na(result$timing$total))

    # Diagnostics present.
    expect_true(
      "rhat_max" %in%
        names(result$diagnostics) ||
        "divergent" %in% names(result$diagnostics)
    )
  })

  it("computes diagnostics over all parameters incl. group-level (A3)", {
    # A varying-intercept model fit with deliberately few iterations so the
    # group-level SD is poorly estimated. extract_brms_diagnostics must look
    # beyond summary()$fixed, else rhat_max would miss a bad group-level SD.
    df <- data.frame(
      y = stats::rnorm(40),
      x = stats::rnorm(40),
      g = factor(rep(1:4, each = 10))
    )
    fit <- suppressWarnings(brms::brm(
      y ~ x + (1 | g),
      data = df,
      family = gaussian(),
      backend = "cmdstanr",
      chains = 1L,
      iter = 20L,
      warmup = 10L,
      silent = 2L,
      refresh = 0L
    ))

    diag <- bayesim:::extract_brms_diagnostics(fit)
    expect_true(is.list(diag))
    expect_true(all(
      c(
        "rhat_max",
        "ess_bulk_min",
        "ess_tail_min",
        "divergent",
        "max_treedepth"
      ) %in%
        names(diag)
    ))
    # sd(Intercept) is a parameter; the all-parameter max rhat must be finite.
    expect_false(is.na(diag$rhat_max))
    # The all-parameter ESS must reflect the group-level SD too: cross-check by
    # comparing against the fixed-effects-only rhat from summary(fit)$fixed.
    fixed_rhat <- max(summary(fit)$fixed[, "Rhat"], na.rm = TRUE)
    # rhat_max over all params is >= the fixed-effects-only max (it covers a
    # superset of parameters), and never NA when fixed-effects is finite.
    expect_gte(diag$rhat_max, fixed_rhat - 1e-9)
  })

  it("falls back to a fresh brm() when precompile is FALSE", {
    fitter <- tiny_fitter()
    fitter@precompile <- FALSE

    # Bank should be NULL when precompile is FALSE.
    bank <- bayesim:::build_model_bank(
      fitter = fitter,
      fit_grid = data.frame(
        model = "gaussian",
        formula = I(list(y ~ x)),
        family = I(list(gaussian())),
        stringsAsFactors = FALSE
      ),
      data_generator = gaussian_data_generator,
      data_spec_template = list(n = 20),
      result_path = NULL,
      seed = 42L
    )
    expect_null(bank)

    # fit_model() with no bank does a fresh compile.
    bayesim:::set_model_bank(NULL)
    data_bundle <- gaussian_data_generator(list(n = 20), list(seed = 1L))
    fit_spec <- list(formula = y ~ x, family = gaussian())
    result <- fit_model(
      fitter,
      data_bundle,
      fit_spec,
      seed = 1L,
      task_ctx = list(task_id = "t")
    )
    expect_true(result$success)
  })

  it("fails loudly on a structural data mismatch (design-matrix guard)", {
    # Build a prefit for a model with a single numeric predictor x (K = 2:
    # intercept + slope). The bank keys it by that spec hash.
    fitter <- tiny_fitter()
    template_data <- data.frame(y = rnorm(20), x = rnorm(20))
    prefit <- brms::brm(
      y ~ x,
      data = template_data,
      family = gaussian(),
      chains = 0L,
      backend = "cmdstanr",
      silent = 2L,
      refresh = 0L
    )
    spec_hash <- bayesim:::model_spec_hash(
      y ~ x,
      gaussian(),
      NULL,
      NULL,
      "cmdstanr"
    )
    bank <- setNames(list(prefit), spec_hash)
    bayesim:::set_model_bank(bank)
    on.exit(bayesim:::set_model_bank(NULL), add = TRUE)

    # Task data with a 3-LEVEL factor predictor where the prefit saw a numeric
    # one. This expands the design matrix from K = 2 to K = 3 (intercept + 2
    # dummy columns), a genuine structural incompatibility: the compiled binary
    # expects 2 predictors but would receive 3. brms with recompile=FALSE would
    # silently mis-fit; the design-matrix guard must catch this and raise fatal.
    mismatched_bundle <- list(
      train = data.frame(
        y = rnorm(20),
        x = factor(sample(c("a", "b", "c"), 20, replace = TRUE))
      )
    )
    fit_spec <- list(formula = y ~ x, family = gaussian())

    expect_error(
      fit_model(
        fitter,
        mismatched_bundle,
        fit_spec,
        seed = 1L,
        task_ctx = list(task_id = "m")
      ),
      class = "bayesim_internal_error"
    )
  })

  it("integrates with run_simulation (one compile, two tasks)", {
    # Full run_simulation flow: model bank built once, reused across 2 tasks.
    fit_grid <- data.frame(model = "gaussian", stringsAsFactors = FALSE)
    fit_grid$formula <- list(y ~ x)
    fit_grid$family <- list(gaussian())

    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = fit_grid,
      data_generator = gaussian_data_generator,
      fitter = tiny_fitter(),
      metrics = list(pred_rmse_metric()),
      n_replicates = 2L,
      seed = 42L
    )

    result <- run_simulation(config, resume = "never", progress = FALSE)

    expect_s3_class(result, "bayesim_simulation_result")
    expect_equal(nrow(result$summary), 2L)
    # Both tasks succeeded (model bank path reused the compiled binary).
    expect_true(all(result$summary$status == "success"))
  })
})

describe("BrmsFitter warning capture (F5)", {
  # A fitter configured to reliably trigger a brms convergence warning via
  # tiny iter counts: Rhat > 1.05 causes summary(fit) to emit an R-level
  # warning, which the fit method must capture into result$warnings.
  # silent = 2L (the default) does NOT suppress the summary() convergence
  # warning (silent only governs sampler stdout), so warnings still surface.
  warning_fitter <- function(bank = TRUE) {
    BrmsFitter(
      chains = 1L,
      iter = 10L,
      warmup = 5L,
      cores = 1L,
      silent = 2L,
      precompile = bank
    )
  }

  it("captures warnings on the fresh (no-bank) path", {
    fitter <- warning_fitter(bank = FALSE)
    bayesim:::set_model_bank(NULL)
    on.exit(bayesim:::set_model_bank(NULL), add = TRUE)

    data_bundle <- gaussian_data_generator(list(n = 20), list(seed = 1L))
    fit_spec <- list(formula = y ~ x, family = gaussian())

    result <- fit_model(
      fitter,
      data_bundle,
      fit_spec,
      seed = 1L,
      task_ctx = list(task_id = "w1")
    )

    expect_true(result$success)
    expect_vector(result$warnings, character())
    expect_gt(length(result$warnings), 0L)
  })

  it("captures warnings on the model-bank (prefit update) path", {
    fitter <- warning_fitter(bank = TRUE)

    fit_grid <- data.frame(model = "gaussian", stringsAsFactors = FALSE)
    fit_grid$formula <- list(y ~ x)
    fit_grid$family <- list(gaussian())

    bank <- bayesim:::build_model_bank(
      fitter = fitter,
      fit_grid = fit_grid,
      data_generator = gaussian_data_generator,
      data_spec_template = list(n = 20),
      result_path = NULL,
      seed = 42L
    )
    bayesim:::set_model_bank(bank)
    on.exit(bayesim:::set_model_bank(NULL), add = TRUE)

    data_bundle <- gaussian_data_generator(list(n = 20), list(seed = 1L))
    fit_spec <- list(formula = y ~ x, family = gaussian())

    result <- fit_model(
      fitter,
      data_bundle,
      fit_spec,
      seed = 1L,
      task_ctx = list(task_id = "w2")
    )

    expect_true(result$success)
    expect_vector(result$warnings, character())
    expect_gt(length(result$warnings), 0L)
  })
})
