# Contributing to bayesim

Thanks for helping improve bayesim. Bug reports and focused pull requests are
welcome.

## Development setup

Install the package dependencies, including development tools, then run:

```r
devtools::document()
devtools::test()
devtools::check()
```

## Test tiers

Tests are split into three explicit tiers, declared in
`tools/ci/test-tiers.R`. Every `tests/testthat/test-*.R` file belongs to
exactly one tier; a `stopifnot` assertion in that script fails CI when a new
file has not been classified, so nothing can silently drop out of every tier.

- **fast** — the fail-fast PR gate. Stan-free analytic workflows (the golden
  complete-study and analytic SBC tests, lifecycle parity, checkpoint/resume,
  external extension contracts, metric schema conformance, statistical
  formula parity) plus focused invariants. Runs in well under 90 seconds.
  CI: `test-fast` on every push and pull request.
- **core** — fast plus the rest of the analytic suite: engine internals,
  analysis layer, model grid, parquet summaries, parallel transport, and
  report rendering. This is what `R CMD check` and a plain `devtools::test()`
  run. CI: `test-coverage`.
- **backend** — requires brms, cmdstanr, and a compiled CmdStan
  (`brms-fitter`, `cmdstan-fitter`, `generators`, `loo-metrics-parity`,
  `metrics-brms`, `sbc-acceptance`). These tests skip unless
  `BAYESIM_TEST_TIER=backend` is set. CI: `test-nightly` and the
  `backend-acceptance` job of `release-check` (one cached Linux runner).

Run a tier locally either through the helper (it applies the file filter and
exports the tier env var for you):

```r
source("tools/ci/test-tiers.R")
run_bayesim_test_tier("fast")
run_bayesim_test_tier("core")
run_bayesim_test_tier("backend")
```

or directly with `devtools::test()` and the tier env var, which is how the
backend-gated skips are controlled:

```sh
BAYESIM_TEST_TIER=fast    Rscript -e 'devtools::test(filter = "^sbc-analytic$")'
BAYESIM_TEST_TIER=backend Rscript -e 'devtools::test()'
```

An unset `BAYESIM_TEST_TIER` means core: every analytic test runs, and
backend-tier tests skip. Run the backend tier locally when you change brms,
cmdstanr, model-bank, transport, or SBC code.

The Stan-backed SBC acceptance suite additionally requires
`BAYESIM_RUN_ACCEPTANCE=true`; without the backend tier set as well, every
test in it skips:

```sh
BAYESIM_TEST_TIER=backend BAYESIM_RUN_ACCEPTANCE=true \
  Rscript -e 'devtools::test(filter = "sbc-acceptance")'
```

## Benchmarks

`inst/benchmarks/checkpoint-scaling.R` measures the controller-side
checkpoint write/read path: it pushes outcomes built with `new_task_result()`
through the sharded RunStore in replicate batches and reports wall time,
outcome shard count, persisted outcome count, checkpoint size, and peak
memory. Run it (from the package root) when you touch `R/checkpoint.R`,
`R/run-store.R`, `R/resume.R`, or the batching logic in `R/simulate.R`, and
compare before/after output; it is the quickest signal for accidental
full-history rewrites or dropped outcomes:

```sh
Rscript inst/benchmarks/checkpoint-scaling.R 1000 10000
BAYESIM_BENCH_BATCH=50 Rscript inst/benchmarks/checkpoint-scaling.R 1000
```

The script verifies its own integrity (persisted outcomes must equal the task
count) and warns if the checkpoint path starts dropping outcomes.

Format R code with `air format .` and run `jarl check .` before submitting a
pull request. Generated documentation (`man/`, `NAMESPACE`, and `README.md`)
should be updated whenever their roxygen or `README.Rmd` sources change.

## Pull requests

- Add a behavioral regression test for bug fixes.
- Keep public behavior and documentation synchronized.
- Explain statistical assumptions and cite the reference implementation or
  paper when adding a statistical measure.
- Avoid committing compiled Stan binaries or local editor/agent state.
- Classify every new `test-*.R` file into a tier in `tools/ci/test-tiers.R`
  (CI fails otherwise).

## Known check warnings

`R CMD check` reports one codoc WARNING for the metric class constructors
(`BiasMetric`, `CoverageMetric`, ...). S7 injects property defaults such as
`schema = list(...)` into the auto-generated constructor formals as
*evaluated values*, while the Rd usage section documents them as
*expressions*; codoc compares the two with `as.character()`, which renders
a list value and an equivalent expression differently. No Rd-level fix
exists — the warning disappears only if S7 starts storing expression
defaults or codoc special-cases value formals. Do not try to silence it by
hand-editing man pages; regenerate docs with `devtools::document()`.
