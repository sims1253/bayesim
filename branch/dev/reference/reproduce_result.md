# This method will reproduce the exact dataset and fit corresponding to the supplied result dataframe row.

The code in this function is written so that all seeds are set at the
right time and all following code after setting the seed replicates
exactly as during the simulation.

## Usage

``` r
reproduce_result(result, data_generator_fn = NULL)
```

## Arguments

- result:

  A data.frame row containing simulation result metadata

- data_generator_fn:

  A function with signature `(config, seed) -> list` that generates
  data. The returned list must contain: `dataset`, `sampling_loops`,
  `bad_samples`, `testing_data`, and `true_parameters`.

## Value

A list containing the fitted model, dataset, and metadata

## Examples

``` r
if (FALSE) { # \dontrun{
# reproduce_result() requires a custom data_generator_fn
# See package vignettes for complete usage examples
} # }
```
