# Explicit CI test tiers. Keep the fast list small: each entry protects a
# complete workflow or a focused invariant that is hard to diagnose end to end.
bayesim_test_tier <- function(tier = c("fast", "backend")) {
  tier <- match.arg(tier)
  switch(
    tier,
    fast = c(
      "adaptive-stop",
      "core-invariants",
      "custom-fitter-public",
      "fitter-orientation",
      "lifecycle-integration",
      "linear-regression-fitter",
      "performance-measures",
      "rng-session",
      "runtime-ux"
    ),
    backend = c(
      "brms-fitter",
      "cmdstan-fitter",
      "generators",
      "loo-metrics-parity",
      "metrics-brms",
      "rstar-metric",
      "sbc-acceptance"
    )
  )
}

bayesim_test_filter <- function(tier = c("fast", "backend")) {
  paste0("^(", paste(bayesim_test_tier(tier), collapse = "|"), ")$")
}

run_bayesim_test_tier <- function(
  tier = c("fast", "backend"),
  reporter = "summary"
) {
  tier <- match.arg(tier)
  withr::local_envvar(BAYESIM_TEST_TIER = tier)
  testthat::test_local(
    filter = bayesim_test_filter(tier),
    reporter = reporter
  )
}
