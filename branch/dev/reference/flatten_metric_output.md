# Flatten Metric Output

Flattens a metric output list into a single-level named list with
prefixed column names. This is used by the simulation engine to convert
metric outputs into flat columns for results collection.

## Usage

``` r
flatten_metric_output(output, metric_name)
```

## Arguments

- output:

  The metric output list to flatten.

- metric_name:

  Character string identifying the metric, used as prefix.

## Value

A named list with flattened names. Scalar values get names like
`<metric_name>__<field>`. Named numeric vectors get expanded to
`<metric_name>__<field>__<subname>`.

## Flattening Rules

- Scalar values become `<metric_name>__<field>`

- Named numeric vectors expand to `<metric_name>__<field>__<subname>`
  for each element

- The double underscore `__` is used as a separator consistently

## Examples

``` r
# Scalar values
flatten_metric_output(list(rmse = 0.5, n_obs = 100L), "my_metric")
#> $my_metric__rmse
#> [1] 0.5
#> 
#> $my_metric__n_obs
#> [1] 100
#> 
# Returns: list(my_metric__rmse = 0.5, my_metric__n_obs = 100L)

# Named numeric vectors are expanded
flatten_metric_output(
  list(params = c(alpha = 0.1, beta = 0.2, gamma = 0.3)),
  "estimates"
)
#> $estimates__params__alpha
#> [1] 0.1
#> 
#> $estimates__params__beta
#> [1] 0.2
#> 
#> $estimates__params__gamma
#> [1] 0.3
#> 
# Returns: list(
#   estimates__params__alpha = 0.1,
#   estimates__params__beta = 0.2,
#   estimates__params__gamma = 0.3
# )
```
