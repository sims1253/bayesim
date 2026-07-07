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

  An S7 SimulationConfig object.

## Value

A character string containing the SHA256 hash of the configuration.

## Details

The fingerprint excludes runtime policy settings (B4):

- `result_path`: Output location doesn't affect simulation identity

- `checkpoint_every` / `checkpoint_format`: checkpoint cadence/format is
  runtime optimization

- `retain`: retention policy must not invalidate resume

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
