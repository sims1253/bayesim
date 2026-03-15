# Access metrics with string identifiers

This function is mostly a helper function that maps identifier strings
to metrics for convenient use in some contexts.

## Usage

``` r
metric_lookup(
  metric,
  fit = NULL,
  draws = NULL,
  testing_data = NULL,
  vars_of_interest = NULL,
  references = NULL,
  threshold = 0.7,
  psis_object = NULL,
  ppred = NULL,
  quantiles = NULL,
  data_gen_output = NULL,
  fit_conf = NULL,
  ...
)
```

## Arguments

- metric:

  A string that identifies a supported metric

- fit:

  A brms fit object

- draws:

  Posterior draws object

- testing_data:

  Data used for testing/out-of-sample evaluation

- vars_of_interest:

  Variables to compute metrics for

- references:

  Reference values for computing distance metrics

- threshold:

  Threshold for diagnostic checks (default: 0.7)

- psis_object:

  PSIS-LOO object for importance sampling

- ppred:

  Posterior predictive samples matrix

- quantiles:

  Quantiles to compute for summaries

- data_gen_output:

  Output from data generation process

- fit_conf:

  Fit configuration list

- ...:

  Additional arguments passed to metric functions

## Value

The function corresponding to the identifier string.
