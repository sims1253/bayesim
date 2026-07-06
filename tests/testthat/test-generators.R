skip_on_cran()

skip_if_not(requireNamespace("cmdstanr", quietly = TRUE))
skip_if_not(requireNamespace("brms", quietly = TRUE))
skip_if_not(requireNamespace("posterior", quietly = TRUE))

# A simple predictor generator that consumes the ambient RNG state.
gaussian_predictors <- function(data_spec, task_ctx) {
  n <- as.integer(data_spec$n %||% 20L)
  data.frame(x = stats::rnorm(n))
}

describe("fixed_truth_generator", {
  it("pins true_params and delegates data generation", {
    gen <- fixed_truth_generator(
      truth = c(beta = 1, sigma = 1),
      draw_data = function(data_spec, task_ctx) {
        n <- as.integer(data_spec$n %||% 20L)
        x <- stats::rnorm(n)
        y <- x + stats::rnorm(n)
        list(train = data.frame(y = y, x = x), response = "y")
      }
    )
    bundle <- gen(list(n = 20), list(task_id = "t1", rep_idx = 1L))
    expect_equal(bundle$true_params, c(beta = 1, sigma = 1))
    expect_equal(bundle$vars_of_interest, c("beta", "sigma"))
    expect_s3_class(bundle$train, "data.frame")
    expect_equal(bundle$meta$generator, "fixed_truth")
  })

  it("rejects non-numeric or unnamed truth", {
    expect_error(
      fixed_truth_generator(c(1, 2), \(d, t) list()),
      class = "bayesim_config_error"
    )
    expect_error(
      fixed_truth_generator(c(a = "x"), \(d, t) list()),
      class = "bayesim_config_error"
    )
  })

  it("rejects draw_data without the right signature", {
    expect_error(
      fixed_truth_generator(c(a = 1), \(x) list()),
      class = "bayesim_config_error"
    )
  })
})

describe("prior_predictive_generator", {
  it("draws theta from the prior and simulates y|theta deterministically", {
    # Compile a sample_prior="only" model once for the block.
    prior_fit <- brms::brm(
      y ~ x,
      data = data.frame(y = rnorm(20), x = rnorm(20)),
      family = gaussian(),
      prior = brms::prior(normal(0, 5), class = "b"),
      sample_prior = "only",
      backend = "cmdstanr",
      chains = 1L,
      iter = 50L,
      warmup = 25L,
      silent = 2L,
      refresh = 0L
    )

    gen <- prior_predictive_generator(
      prior_fit = prior_fit,
      predictor_generator = gaussian_predictors,
      vars_of_interest = "x"
    )

    # Same rep_idx -> identical true_params (the fixed prior draw theta).
    b1 <- gen(list(n = 20), list(task_id = "a", rep_idx = 1L))
    b2 <- gen(list(n = 20), list(task_id = "b", rep_idx = 1L))
    expect_equal(b1$true_params, b2$true_params)

    # The simulated y depends on the ambient RNG (for predictors), so it only
    # matches when the RNG state is restored; the draw-indexed theta is the
    # deterministic signal. Verify truth_draw_id is fixed by rep_idx.
    expect_equal(b1$meta$truth_draw_id, b2$meta$truth_draw_id)

    # Different rep_idx -> different theta (when more than 1 prior draw).
    if (posterior::ndraws(prior_fit) > 1L) {
      b3 <- gen(list(n = 20), list(task_id = "c", rep_idx = 2L))
      expect_false(identical(b1$true_params, b3$true_params))
    }

    expect_equal(b1$meta$generator, "prior_predictive")
    expect_true(!is.null(b1$meta$truth_draw_id))
    expect_equal(b1$meta$truth_draw_id, 1L)
  })
})

describe("ifs_generator", {
  it("draws theta from a preconditioning fit and forward-simulates", {
    # A real preconditioning fit (posterior, not prior-only).
    prefit <- brms::brm(
      y ~ x,
      data = data.frame(y = rnorm(20), x = rnorm(20)),
      family = gaussian(),
      backend = "cmdstanr",
      chains = 1L,
      iter = 50L,
      warmup = 25L,
      silent = 2L,
      refresh = 0L
    )

    gen <- ifs_generator(
      prefit = prefit,
      predictor_generator = gaussian_predictors,
      vars_of_interest = "x"
    )

    b1 <- gen(list(n = 20), list(task_id = "a", rep_idx = 1L))
    b2 <- gen(list(n = 20), list(task_id = "b", rep_idx = 1L))
    expect_equal(b1$true_params, b2$true_params)

    expect_equal(b1$meta$generator, "ifs")
    expect_equal(b1$meta$truth_draw_id, 1L)
    expect_s3_class(b1$train, "data.frame")
  })

  it("applies truncation bounds when requested", {
    prefit <- brms::brm(
      y ~ x,
      data = data.frame(y = rnorm(20), x = rnorm(20)),
      family = gaussian(),
      backend = "cmdstanr",
      chains = 1L,
      iter = 50L,
      warmup = 25L,
      silent = 2L,
      refresh = 0L
    )
    gen <- ifs_generator(
      prefit = prefit,
      predictor_generator = gaussian_predictors,
      vars_of_interest = "x",
      lower_bound = -1e6, # effectively no truncation, but exercises the path
      upper_bound = 1e6,
      truncate = TRUE
    )
    b <- gen(list(n = 20), list(task_id = "t", rep_idx = 1L))
    expect_true(all(b$train$y >= -1e6))
    expect_true(all(b$train$y <= 1e6))
  })
})

describe("brms_response_sequence internals", {
  it("returns list(y) for a univariate model (F1 regression)", {
    prefit <- brms::brm(
      y ~ x,
      data = data.frame(y = rnorm(20), x = rnorm(20)),
      family = gaussian(),
      backend = "cmdstanr",
      chains = 1L,
      iter = 50L,
      warmup = 25L,
      silent = 2L,
      refresh = 0L
    )
    seq <- bayesim:::brms_response_sequence(prefit)
    expect_type(seq, "list")
    expect_length(seq, 1L)
    # One depth group containing the response name "y" — NOT a dependency list.
    expect_equal(seq, list("y"))
  })

  it("orders responses by dependency for a multivariate model", {
    # b depends on a (b's formula includes a); a depends only on x.
    prefit <- brms::brm(
      brms::bf(a ~ x) + brms::bf(b ~ a),
      data = data.frame(a = rnorm(20), b = rnorm(20), x = rnorm(20)),
      family = gaussian(),
      backend = "cmdstanr",
      chains = 1L,
      iter = 50L,
      warmup = 25L,
      silent = 2L,
      refresh = 0L
    )
    seq <- bayesim:::brms_response_sequence(prefit)
    expect_type(seq, "list")
    expect_gte(length(seq), 1L)
    # The union of all depth groups covers both responses.
    all_resp <- unique(unlist(seq))
    expect_setequal(all_resp, c("a", "b"))
    # a must be simulated before b (a appears in an earlier-or-equal depth group
    # index than b). For this model they land in separate groups.
    depth_of <- function(name) {
      min(which(vapply(seq, function(g) name %in% g, logical(1))))
    }
    expect_equal(depth_of("a"), 1L)
    expect_equal(depth_of("b"), 2L)
  })

  it("errors on a cyclic response dependency (nodes_by_depth)", {
    # A 2x2 adjacency matrix where each node points to the other (a cycle):
    # rows have nonzero sums, so no zero-indegree node exists.
    cyclic <- matrix(
      c(0, 1, 1, 0),
      nrow = 2,
      ncol = 2,
      dimnames = list(c("a", "b"), c("a", "b"))
    )
    expect_error(
      bayesim:::nodes_by_depth(cyclic),
      class = "bayesim_config_error"
    )
  })
})

describe("ifs_generator F1 — simulated response", {
  it("simulates a y column that differs from the pilot data (predictor_generator supplied)", {
    prefit <- brms::brm(
      y ~ x,
      data = data.frame(y = rnorm(20), x = rnorm(20)),
      family = gaussian(),
      backend = "cmdstanr",
      chains = 1L,
      iter = 50L,
      warmup = 25L,
      silent = 2L,
      refresh = 0L
    )
    gen <- ifs_generator(
      prefit = prefit,
      predictor_generator = gaussian_predictors,
      vars_of_interest = "x"
    )
    bundle <- gen(list(n = 20), task_ctx = list(task_id = "t1", rep_idx = 1L))
    expect_true("y" %in% names(bundle$train))
    expect_false(all(is.na(bundle$train$y)))
    # The simulated y must NOT be literally the pilot data's y (it is freshly
    # drawn from the posterior predictive).
    expect_false(identical(bundle$train$y, prefit$data$y))
  })

  it("simulates a y column even when predictor_generator = NULL (F1 regression)", {
    prefit <- brms::brm(
      y ~ x,
      data = data.frame(y = rnorm(20), x = rnorm(20)),
      family = gaussian(),
      backend = "cmdstanr",
      chains = 1L,
      iter = 50L,
      warmup = 25L,
      silent = 2L,
      refresh = 0L
    )
    gen <- ifs_generator(
      prefit = prefit,
      predictor_generator = NULL,
      vars_of_interest = "x"
    )
    bundle <- gen(list(n = 20), task_ctx = list(task_id = "t1", rep_idx = 1L))
    # Predictor columns come from the pilot, but the response MUST be freshly
    # simulated (previously this returned the pilot frame unchanged).
    expect_true("y" %in% names(bundle$train))
    expect_false(all(is.na(bundle$train$y)))
    expect_false(identical(bundle$train$y, prefit$data$y))
  })

  it("produces a valid bundle for rep_idx = 2 and a value > n_draws (modulo wrap)", {
    prefit <- brms::brm(
      y ~ x,
      data = data.frame(y = rnorm(20), x = rnorm(20)),
      family = gaussian(),
      backend = "cmdstanr",
      chains = 1L,
      iter = 50L,
      warmup = 25L,
      silent = 2L,
      refresh = 0L
    )
    n_draws <- posterior::ndraws(prefit)
    gen <- ifs_generator(
      prefit = prefit,
      predictor_generator = gaussian_predictors,
      vars_of_interest = "x"
    )

    b2 <- gen(list(n = 20), task_ctx = list(task_id = "t2", rep_idx = 2L))
    expect_equal(b2$meta$truth_draw_id, 2L)
    expect_true("y" %in% names(b2$train))
    expect_false(all(is.na(b2$train$y)))

    # rep_idx beyond n_draws must wrap via modulo without indexing errors.
    rep_idx_wrap <- as.integer(n_draws) + 5L
    expected_id <- ((rep_idx_wrap - 1L) %% n_draws) + 1L
    bw <- gen(
      list(n = 20),
      task_ctx = list(task_id = "tw", rep_idx = rep_idx_wrap)
    )
    expect_equal(bw$meta$truth_draw_id, expected_id)
    expect_false(all(is.na(bw$train$y)))
  })

  it("simulates BOTH response columns for a multivariate model", {
    prefit <- brms::brm(
      brms::bf(a ~ x) + brms::bf(b ~ a),
      data = data.frame(a = rnorm(20), b = rnorm(20), x = rnorm(20)),
      family = gaussian(),
      backend = "cmdstanr",
      chains = 1L,
      iter = 50L,
      warmup = 25L,
      silent = 2L,
      refresh = 0L
    )
    gen <- ifs_generator(
      prefit = prefit,
      predictor_generator = gaussian_predictors
    )
    bundle <- gen(list(n = 20), task_ctx = list(task_id = "mv1", rep_idx = 1L))
    expect_true("a" %in% names(bundle$train))
    expect_true("b" %in% names(bundle$train))
    expect_false(all(is.na(bundle$train$a)))
    expect_false(all(is.na(bundle$train$b)))
  })
})

describe("prior_predictive_generator F1 — simulated response", {
  it("fills the response column and varies it across draws", {
    prior_fit <- brms::brm(
      y ~ x,
      data = data.frame(y = rnorm(20), x = rnorm(20)),
      family = gaussian(),
      prior = brms::prior(normal(0, 5), class = "b"),
      sample_prior = "only",
      backend = "cmdstanr",
      chains = 1L,
      iter = 50L,
      warmup = 25L,
      silent = 2L,
      refresh = 0L
    )
    gen <- prior_predictive_generator(
      prior_fit = prior_fit,
      predictor_generator = gaussian_predictors,
      vars_of_interest = "x"
    )
    b1 <- gen(list(n = 20), list(task_id = "a", rep_idx = 1L))
    # Response column present and non-NA.
    expect_true("y" %in% names(b1$train))
    expect_false(all(is.na(b1$train$y)))
    # Different rep_idx -> different simulated response (varies with the draw).
    if (posterior::ndraws(prior_fit) > 1L) {
      b2 <- gen(list(n = 20), list(task_id = "b", rep_idx = 2L))
      expect_false(identical(b1$train$y, b2$train$y))
    }
  })
})
