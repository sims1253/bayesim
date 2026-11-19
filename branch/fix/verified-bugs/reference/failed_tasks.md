# Extract failed tasks from a simulation result

Returns the failed-task errors joined with grid columns, as a tibble —
convenient for diagnosing what went wrong after a run with failures.
Also used by run_simulation to print a compact failure summary.

## Usage

``` r
failed_tasks(result)
```

## Arguments

- result:

  A `bayesim_simulation_result`.

## Value

A tibble with grid columns plus `error_class`, `error_message`.

## See also

[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
failed_tasks(result)
} # }
```
