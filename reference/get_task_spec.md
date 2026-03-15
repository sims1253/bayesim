# Get Task Specification from Grid

Extracts the complete specification for a single task, including
data_spec, fit_spec, task context, and RNG seed.

## Usage

``` r
get_task_spec(task_grid, task_id, config)
```

## Arguments

- task_grid:

  A task grid tibble created by
  [`create_task_grid()`](https://sims1253.github.io/bayesim/reference/create_task_grid.md).

- task_id:

  The task ID to look up (format: "dXXX_fXXX_rXXXXX").

- config:

  The SimulationConfig object used to create the task grid.

## Value

A named list with elements:

- `task_id`: Character task identifier

- `data_idx`: Integer data grid row index

- `fit_idx`: Integer fit grid row index

- `rep_idx`: Integer replicate number

- `data_spec`: Named list of data specification parameters

- `fit_spec`: Named list of fit specification parameters

- `task_ctx`: Named list with task_id, data_idx, fit_idx, rep_idx

- `rng_seed`: Integer vector RNG state for this task

## Examples

``` r
if (FALSE) { # \dontrun{
task_grid <- create_task_grid(config)
spec <- get_task_spec(task_grid, "d001_f001_r00001", config)
spec$data_spec  # Data parameters for this task
spec$fit_spec   # Fit parameters for this task
} # }
```
