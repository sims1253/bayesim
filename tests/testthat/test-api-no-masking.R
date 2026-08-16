# B1: the renamed generics must not mask common foreign generics.
# bayesim no longer exports fit / compute / log_lik / diagnostics, so loading it
# alongside brms/dplyr must leave those names resolving to the foreign package.
# This file is fast-tier analytic; the brms dispatch companion lives in
# test-brms-fitter.R (backend tier) so it actually executes.

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
})
