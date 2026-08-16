# Resume a simulation from an existing result directory

Calls
[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md)
with `resume = "must"`.

## Usage

``` r
resume_simulation(result_path, config = NULL, progress = TRUE, verbose = TRUE)
```

## Arguments

- result_path:

  Character path to an existing result directory

- config:

  Optional SimulationConfig object. When omitted, bayesim rebuilds the
  config from the run manifest. That works only when every executable
  component (data generator, fitter, metrics) is a namespaced package
  function or class. The run manifest cannot rehydrate script-defined
  closures, such as the return value of
  [`fixed_truth_generator()`](https://sims1253.github.io/bayesim/reference/fixed_truth_generator.md)
  or any generator defined in your script: R can serialize closures, but
  configless resume intentionally does not restore arbitrary executable
  closures. Those runs require the original `config`.

- progress:

  Logical; if TRUE, show progress bar

- verbose:

  Logical; if TRUE, print lifecycle messages independently of the
  progress bar.

## Value

A bayesim_simulation_result S3 object
