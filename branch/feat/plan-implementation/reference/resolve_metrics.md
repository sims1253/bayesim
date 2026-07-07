# Resolve Metrics to List of Metric Objects

Internal helper to convert character metric names or Metric objects into
a standardized list of Metric objects.

## Usage

``` r
resolve_metrics(metrics)
```

## Arguments

- metrics:

  Character vector of metric names, list of Metric objects, or NULL.

## Value

A list of Metric objects, or NULL if input was NULL.
