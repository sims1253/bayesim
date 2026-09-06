# Get Total Number of Tasks in Configuration

Calculates the total number of simulation tasks based on the
configuration. Each task is one (data_spec, fit_spec, replicate)
combination.

## Usage

``` r
get_total_tasks(config)
```

## Arguments

- config:

  An S7 SimulationConfig object.

## Value

Integer. Total number of tasks.

## Examples

``` r
if (FALSE) { # \dontrun{
config <- simulation_config(
  data_grid = data.frame(n = c(100, 500)),  # 2 rows
  fit_grid = data.frame(model = c("A", "B")),  # 2 rows
  n_replicates = 100L,
  ...
)
get_total_tasks(config)  # Returns 400 (2 * 2 * 100)
} # }
```
