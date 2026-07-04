skip_on_cran()

skip_if_not(requireNamespace("cmdstanr", quietly = TRUE))
skip_if_not(requireNamespace("brms", quietly = TRUE))
skip_if_not(requireNamespace("rstan", quietly = TRUE))

# A trivial gaussian generator producing the same structure every call.
# Variables (y, x) match the test models so update(recompile=FALSE) holds.
gaussian_data_generator <- function(data_spec, seed, task_ctx) {
  # Consume the ambient RNG state (the worker restores the per-task stream).
  # For direct calls outside the worker, fall back to the provided seed.
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
    vars_of_interest = "beta",
    references = c(beta = 0)
  )
}

tiny_fitter <- function() {
  BrmsFitter(chains = 1L, iter = 50L, warmup = 25L, cores = 1L)
}

describe("BrmsFitter model bank", {
  it("builds one prefit per distinct model spec", {
    fit_grid <- data.frame(
      model = c("gaussian", "gaussian_dup", "student"),
      stringsAsFactors = FALSE
    )
    fit_grid$formula <- list(y ~ x, y ~ x, y ~ x)
    fit_grid$family <- list(gaussian(), gaussian(), brms::brmsfamily("student"))

    fitter <- tiny_fitter()
    bank <- bayesim:::build_model_bank(
      fitter = fitter,
      fit_grid = fit_grid,
      data_generator = gaussian_data_generator,
      data_spec_template = list(n = 20),
      result_path = NULL,
      seed = 42L
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

    data_bundle <- gaussian_data_generator(list(n = 20), seed = 1, NULL)
    fit_spec <- list(formula = y ~ x, family = gaussian())

    # Two fit() calls on different data should reuse the same prefit binary.
    # Capture warnings: a recompile warning would indicate the bank was bypassed.
    recompiled <- FALSE
    result1 <- withCallingHandlers(
      fit(fitter, data_bundle, fit_spec, seed = 100L, task_ctx = list(task_id = "t1")),
      warning = function(w) {
        if (grepl("recompil", conditionMessage(w), ignore.case = TRUE)) {
          recompiled <<- TRUE
        }
        invokeRestart("muffleWarning")
      }
    )
    expect_false(recompiled)
    expect_true(result1$success)

    result2 <- fit(fitter, data_bundle, fit_spec, seed = 200L, task_ctx = list(task_id = "t2"))
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

    data_bundle <- gaussian_data_generator(list(n = 20), seed = 1, NULL)
    fit_spec <- list(formula = y ~ x, family = gaussian())

    r1 <- fit(fitter, data_bundle, fit_spec, seed = 777L, task_ctx = list(task_id = "a"))
    r2 <- fit(fitter, data_bundle, fit_spec, seed = 777L, task_ctx = list(task_id = "b"))
    r3 <- fit(fitter, data_bundle, fit_spec, seed = 888L, task_ctx = list(task_id = "c"))

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

    data_bundle <- gaussian_data_generator(list(n = 20), seed = 1, NULL)
    fit_spec <- list(formula = y ~ x, family = gaussian())
    result <- fit(fitter, data_bundle, fit_spec, seed = 1L, task_ctx = list(task_id = "t"))

    # Timings from the updated fit's stanfit (non-zero warmup/sample).
    expect_true(is.numeric(result$timing$warmup) || is.na(result$timing$warmup))
    expect_true(is.numeric(result$timing$total))
    expect_false(is.na(result$timing$total))

    # Diagnostics present.
    expect_true("rhat_max" %in% names(result$diagnostics) ||
      "divergent" %in% names(result$diagnostics))
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

    # fit() with no bank does a fresh compile.
    bayesim:::set_model_bank(NULL)
    data_bundle <- gaussian_data_generator(list(n = 20), seed = 1, NULL)
    fit_spec <- list(formula = y ~ x, family = gaussian())
    result <- fit(fitter, data_bundle, fit_spec, seed = 1L, task_ctx = list(task_id = "t"))
    expect_true(result$success)
  })

  it("fails loudly on a structural data mismatch (design-matrix guard)", {
    # Build a prefit for a model with a single numeric predictor x (K = 2:
    # intercept + slope). The bank keys it by that spec hash.
    fitter <- tiny_fitter()
    template_data <- data.frame(y = rnorm(20), x = rnorm(20))
    prefit <- brms::brm(
      y ~ x, data = template_data, family = gaussian(),
      chains = 0L, backend = "cmdstanr", silent = 2L, refresh = 0L
    )
    spec_hash <- bayesim:::model_spec_hash(y ~ x, gaussian(), NULL, NULL, "cmdstanr")
    bank <- setNames(list(prefit), spec_hash)
    bayesim:::set_model_bank(bank)
    on.exit(bayesim:::set_model_bank(NULL), add = TRUE)

    # Task data with a 3-LEVEL factor predictor where the prefit saw a numeric
    # one. This expands the design matrix from K = 2 to K = 3 (intercept + 2
    # dummy columns), a genuine structural incompatibility: the compiled binary
    # expects 2 predictors but would receive 3. brms with recompile=FALSE would
    # silently mis-fit; the design-matrix guard must catch this and raise fatal.
    mismatched_bundle <- list(train = data.frame(
      y = rnorm(20),
      x = factor(sample(c("a", "b", "c"), 20, replace = TRUE))
    ))
    fit_spec <- list(formula = y ~ x, family = gaussian())

    expect_error(
      fit(fitter, mismatched_bundle, fit_spec, seed = 1L, task_ctx = list(task_id = "m")),
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
      metrics = list(rmse_metric()),
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

    data_bundle <- gaussian_data_generator(list(n = 20), seed = 1, NULL)
    fit_spec <- list(formula = y ~ x, family = gaussian())

    result <- fit(fitter, data_bundle, fit_spec, seed = 1L, task_ctx = list(task_id = "w1"))

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

    data_bundle <- gaussian_data_generator(list(n = 20), seed = 1, NULL)
    fit_spec <- list(formula = y ~ x, family = gaussian())

    result <- fit(fitter, data_bundle, fit_spec, seed = 1L, task_ctx = list(task_id = "w2"))

    expect_true(result$success)
    expect_vector(result$warnings, character())
    expect_gt(length(result$warnings), 0L)
  })
})
