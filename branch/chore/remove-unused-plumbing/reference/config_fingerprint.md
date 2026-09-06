# Stable simulation-design fingerprint

Returns the stable hash bayesim uses to validate checkpoint
compatibility. The fingerprint covers the study design (grids,
generator, fitter, metrics, replicate count, and seed) and intentionally
excludes runtime policy such as output paths, retention, checkpoint
cadence, and adaptive stopping.

## Usage

``` r
config_fingerprint(config)
```

## Arguments

- config:

  A simulation configuration created by
  [`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md).

## Value

A scalar SHA-256 character string.

## See also

[`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md),
[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
config_fingerprint(config)
} # }
```
