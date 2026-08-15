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

Tests are split into two explicit tiers, defined in `tools/ci/test-tiers.R`:

```r
source("tools/ci/test-tiers.R")
run_bayesim_test_tier("fast")    # no Stan required; the test-fast CI job
run_bayesim_test_tier("backend") # needs brms, cmdstanr, and CmdStan
```

The `fast` tier covers the analytic engine, contracts, and analysis layer; CI
runs it on every push and pull request. The `backend` tier covers brms,
cmdstanr, model-bank, transport, LOO parity, and SBC acceptance; CI runs it
nightly. Backend-tier tests skip unless `BAYESIM_TEST_TIER=backend` is set, so
a plain `devtools::test()` runs only the no-Stan suites. Run the backend tier
locally when you change brms, cmdstanr, model-bank, transport, or SBC code.

The Stan-backed SBC acceptance suite additionally requires
`BAYESIM_RUN_ACCEPTANCE=true`:

```sh
BAYESIM_RUN_ACCEPTANCE=true Rscript -e 'devtools::test(filter = "sbc-acceptance")'
```

Format R code with `air format .` and run `jarl check .` before submitting a
pull request. Generated documentation (`man/`, `NAMESPACE`, and `README.md`)
should be updated whenever their roxygen or `README.Rmd` sources change.

## Pull requests

- Add a behavioral regression test for bug fixes.
- Keep public behavior and documentation synchronized.
- Explain statistical assumptions and cite the reference implementation or
  paper when adding a statistical measure.
- Avoid committing compiled Stan binaries or local editor/agent state.
