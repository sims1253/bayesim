# Contributing to bayesim

Thanks for helping improve bayesim. Bug reports and focused pull
requests are welcome.

## Development setup

Install the package dependencies, including development tools, then run:

``` r

devtools::document()
devtools::test()
devtools::check()
```

The default test loop skips the slow Stan-backed SBC acceptance suite.
Run it explicitly when changing brms, cmdstanr, model-bank, transport,
or SBC code:

``` sh
BAYESIM_RUN_ACCEPTANCE=true Rscript -e 'devtools::test(filter = "sbc-acceptance")'
```

Format R code with `air format .` and run `jarl check .` before
submitting a pull request. Generated documentation (`man/`, `NAMESPACE`,
and `README.md`) should be updated whenever their roxygen or
`README.Rmd` sources change.

## Pull requests

- Add a behavioral regression test for bug fixes.
- Keep public behavior and documentation synchronized.
- Explain statistical assumptions and cite the reference implementation
  or paper when adding a statistical measure.
- Avoid committing compiled Stan binaries or local editor/agent state.
