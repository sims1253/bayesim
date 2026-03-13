# Bias Metric

Mean bias of predictions.

`bias_metric()` is a constructor function that creates a BiasMetric
instance with appropriate defaults. Use this constructor to work around
S7's property default inheritance issue.

## Usage

``` r
BiasMetric(name = character(0), needs = character(0), required = FALSE)

bias_metric(name = "bias")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "bias".

- needs:

  Character vector of required capabilities. Defaults to "predictions".

- required:

  Logical indicating if metric failure causes task failure. Defaults to
  FALSE.

## Value

A BiasMetric object.

A BiasMetric object.

## Examples

``` r
bias_metric()
#> <bayesim::BiasMetric>
#>  @ name    : chr "bias"
#>  @ needs   : chr "predictions"
#>  @ required: logi FALSE
bias_metric(name = "my_bias")
#> <bayesim::BiasMetric>
#>  @ name    : chr "my_bias"
#>  @ needs   : chr "predictions"
#>  @ required: logi FALSE
```
