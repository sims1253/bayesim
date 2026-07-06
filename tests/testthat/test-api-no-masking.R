# B1: the renamed generics must not mask common foreign generics.
# bayesim no longer exports fit / compute / log_lik / diagnostics, so loading it
# alongside brms/dplyr must leave those names resolving to the foreign package.

skip_on_cran()

describe("B1: renamed generics do not mask foreign packages", {
  it("bayesim does not export fit / compute / log_lik / diagnostics", {
    ns <- asNamespace("bayesim")
    exports <- getNamespaceExports(ns)
    expect_false("fit" %in% exports)
    expect_false("compute" %in% exports)
    expect_false("log_lik" %in% exports)
    expect_false("diagnostics" %in% exports)
    # ...and the new names are exported.
    expect_true("fit_model" %in% exports)
    expect_true("compute_metric" %in% exports)
    expect_true("log_lik_matrix" %in% exports)
    expect_true("fit_diagnostics" %in% exports)
  })

  it("log_lik(brmsfit) still resolves to brms after loading bayesim", {
    skip_if_not(requireNamespace("brms", quietly = TRUE))
    skip_if_not(requireNamespace("cmdstanr", quietly = TRUE))

    fit <- suppressWarnings(brms::brm(
      y ~ x,
      data = data.frame(y = rnorm(20), x = rnorm(20)),
      family = gaussian(),
      backend = "cmdstanr",
      chains = 1L,
      iter = 50L,
      warmup = 25L,
      silent = 2L,
      refresh = 0L
    ))

    # Attach brms then bayesim (the package under test) after it. log_lik must
    # still dispatch to brms because bayesim no longer exports a log_lik generic.
    library(brms)
    on.exit(detach("package:brms"), add = TRUE)
    library(bayesim)
    on.exit(detach("package:bayesim"), add = TRUE)

    ll <- log_lik(fit)
    expect_true(is.matrix(ll))
  })
})
