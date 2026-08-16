# Explicit CI test tiers. Every tests/testthat/test-*.R file must appear in
# exactly one tier below; the stopifnot block at the bottom of this script
# fails loudly when a new file has not been classified, so nothing can silently
# fall out of every tier.
#
# Semantics:
# - fast:    the fail-fast PR gate. Stan-free, quick, and focused on complete
#            analytic workflows plus invariants that are hard to diagnose end
#            to end. CI: test-fast (every push/PR). Keep the total wall time
#            under ~90 seconds.
# - core:    fast + everything else analytic. This is what R CMD check and a
#            plain devtools::test() run (BAYESIM_TEST_TIER unset). Slower or
#            secondary suites (engine internals, parallel transport, report
#            rendering) live here. CI: test-coverage.
# - backend: requires a compiled Stan backend (brms/cmdstanr/CmdStan). Skipped
#            unless BAYESIM_TEST_TIER=backend. CI: test-backend and the
#            release-check backend acceptance job.
bayesim_test_tier_fast <- c(
  "adaptive-stop",
  "api-no-masking",
  "checkpoint",
  "contracts",
  "core-invariants",
  "custom-fitter-public",
  "external-metric-public",
  "fitter-agnostic-generators",
  "fitter-orientation",
  "golden-analytic-study",
  "lifecycle-integration",
  "linear-regression-fitter",
  "memory-bounded",
  "metric-schema-conformance",
  "metrics-extended",
  "performance-measures",
  "power-precision",
  "pred-metrics-no-test",
  "regression",
  "rng-session",
  "runtime-ux",
  "sbc-analytic"
)

bayesim_test_tier_core_only <- c(
  "analysis",
  "engine",
  "model-grid",
  "parquet-summary",
  "report",
  "transport-purrr-mirai"
)

bayesim_test_tier_backend <- c(
  "brms-fitter",
  "cmdstan-fitter",
  "generators",
  "loo-metrics-parity",
  "metrics-brms",
  "sbc-acceptance"
)

bayesim_test_tier <- function(tier = c("fast", "core", "backend")) {
  tier <- match.arg(tier)
  switch(
    tier,
    fast = bayesim_test_tier_fast,
    # The core tier is defined as fast + the rest of the analytic suite: it is
    # the tier a plain package check exercises.
    core = c(bayesim_test_tier_fast, bayesim_test_tier_core_only),
    backend = bayesim_test_tier_backend
  )
}

bayesim_test_filter <- function(tier = c("fast", "core", "backend")) {
  paste0("^(", paste(bayesim_test_tier(tier), collapse = "|"), ")$")
}

run_bayesim_test_tier <- function(
  tier = c("fast", "core", "backend"),
  reporter = "summary"
) {
  tier <- match.arg(tier)
  withr::local_envvar(BAYESIM_TEST_TIER = tier)
  testthat::test_local(
    filter = bayesim_test_filter(tier),
    reporter = reporter
  )
}

# Locate tests/testthat relative to the package root. CI jobs and local use
# source this file from the root; walk upward from the working directory as a
# fallback so it also works from nested paths.
bayesim_test_files_on_disk <- function() {
  root <- getwd()
  while (
    length(root) > 1L &&
      !dir.exists(file.path(root, "tests", "testthat")) &&
      !file.exists(file.path(root, "DESCRIPTION"))
  ) {
    root <- dirname(root)
  }
  tests_dir <- file.path(root, "tests", "testthat")
  if (!dir.exists(tests_dir)) {
    stop("Could not locate tests/testthat from ", getwd())
  }
  sub(
    "\\.R$",
    "",
    sub("^test-", "", list.files(tests_dir, pattern = "^test-.*\\.R$"))
  )
}

# Completeness guard: the three tiers must partition the on-disk test files.
# A new test file that has not been assigned a tier fails this assertion (and
# therefore every CI job that sources this script) until it is classified.
bayesim_assert_tier_coverage <- function() {
  declared <- c(
    fast = bayesim_test_tier_fast,
    core = bayesim_test_tier_core_only,
    backend = bayesim_test_tier_backend
  )
  on_disk <- bayesim_test_files_on_disk()

  unclassified <- setdiff(on_disk, declared)
  missing <- setdiff(declared, on_disk)
  duplicated_names <- declared[duplicated(declared)]

  stopifnot(
    "test files missing from every tier (add them to tools/ci/test-tiers.R)" = length(
      unclassified
    ) ==
      0L,
    "tier entries with no matching test file on disk" = length(missing) == 0L,
    "test files listed in more than one tier" = length(duplicated_names) == 0L
  )
  invisible(TRUE)
}
bayesim_assert_tier_coverage()
