# Build the model bank for a BrmsFitter

For each DISTINCT row of `fit_grid` (deduped by
[`model_spec_hash()`](https://sims1253.github.io/bayesim/reference/model_spec_hash.md)),
compiles a prefit via `brms::brm(chains = 0)` against generator-supplied
template data. Returns a named list of prefit objects keyed by spec
hash.

## Usage

``` r
build_model_bank(
  fitter,
  fit_grid,
  data_generator,
  data_spec_template,
  result_path = NULL,
  seed = NULL
)
```

## Arguments

- fitter:

  A BrmsFitter S7 object.

- fit_grid:

  A data.frame of model fitting specifications. Each row's `formula`,
  `family`, `prior`, `stanvars` columns (list-columns) define a model
  spec.

- data_generator:

  The simulation data generator function.

- data_spec_template:

  A named list (one row of `data_grid` as a list) used to generate
  template data.

- result_path:

  NULL or character path. When non-NULL,
  `options(cmdstanr_write_stan_file_dir)` is set to
  `file.path(result_path, "stan_binaries")` so compiled binaries persist
  and are shared across controller and local daemons.

- seed:

  Integer seed for the simulation (used only for logging).

## Value

A named list of `brmsfit` prefit objects keyed by
[`model_spec_hash()`](https://sims1253.github.io/bayesim/reference/model_spec_hash.md),
or NULL when `precompile` is FALSE.

## Details

When `fitter@precompile` is FALSE, returns NULL immediately (the fitter
falls back to a fresh
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) per
task).

A compile failure is fatal: it raises a
[`bayesim_internal_error()`](https://sims1253.github.io/bayesim/reference/bayesim_internal_error.md)
since the simulation cannot proceed without a compilable model.
