# Compute Configuration Fingerprint

Computes a cryptographic hash of the simulation configuration. The
fingerprint uniquely identifies a simulation configuration for caching
and deduplication purposes.

## Usage

``` r
compute_config_fingerprint(config)
```

## Arguments

- config:

  An S7 SimulationConfig object or a `StudySpec` created by
  [`new_study_spec()`](https://sims1253.github.io/bayesim/reference/new_study_spec.md).

## Value

A character string containing the SHA256 hash of the study design.

## Details

The fingerprint excludes runtime policy settings (B4):

- `result_path`: Output location doesn't affect simulation identity

- `checkpoint_every` / `checkpoint_format`: checkpoint cadence/format is
  runtime optimization

- `retain`: retention is runtime policy; fingerprint exclusion does not
  make every retention change legal on resume (a compatible resume may
  narrow retention; widening is rejected once completed outcomes lack
  the requested artifacts)

- `max_errors`: error tolerance is runtime policy

## Examples

``` r
if (FALSE) { # \dontrun{
config <- simulation_config(
  data_grid = data.frame(n = 100),
  fit_grid = data.frame(model = "baseline"),
  data_generator = my_data_gen,
  n_replicates = 10L,
  seed = 42L
)

fingerprint <- compute_config_fingerprint(config)
# Use fingerprint for caching or deduplication
} # }
```
