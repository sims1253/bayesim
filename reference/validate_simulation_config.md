# Validate Simulation Configuration

Validates that a SimulationConfig object is complete and properly
configured for running a simulation. This includes validating the fitter
interface, all metrics, and the data_generator function signature.

## Usage

``` r
validate_simulation_config(config)
```

## Arguments

- config:

  An S7 SimulationConfig object to validate.

## Value

TRUE if validation passes.

## Details

The configuration is validated for:

- Being a valid SimulationConfig S7 object

- Having a non-NULL fitter that passes
  [`validate_fitter_interface()`](https://sims1253.github.io/bayesim/reference/check_fitter_class.md)

- Having metrics (if present) that pass
  [`validate_metric_interface()`](https://sims1253.github.io/bayesim/reference/validate_metric_interface.md)

- Having a data_generator function with the correct signature

## Errors

Throws a `bayesim_config_error` condition if validation fails.

## See also

[`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md),
[`validate_fitter_interface()`](https://sims1253.github.io/bayesim/reference/check_fitter_class.md),
[`validate_metric_interface()`](https://sims1253.github.io/bayesim/reference/validate_metric_interface.md)

## Examples

``` r
if (FALSE) { # \dontrun{
config <- simulation_config(
  data_grid = data.frame(n = 100),
  fit_grid = data.frame(model = "baseline"),
  data_generator = my_data_gen,
  fitter = my_fitter,
  metrics = list(my_metric),
  n_replicates = 10L,
  seed = 42L
)
validate_simulation_config(config)
} # }
```
