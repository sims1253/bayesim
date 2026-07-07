# Bias Metric

Mean bias of predictions.

Constructor for BiasMetric.

## Usage

``` r
BiasMetric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean"
)

pred_bias_metric(name = "bias")
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
