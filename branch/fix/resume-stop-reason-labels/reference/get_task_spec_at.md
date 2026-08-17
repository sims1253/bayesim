# Get a task specification by precomputed row index

The execution loop resolves a whole batch of task IDs once and calls
this helper by index, avoiding repeated full-grid scans for large
studies.

## Usage

``` r
get_task_spec_at(task_grid, row_idx, config)
```

## Arguments

- task_grid:

  A task grid tibble.

- row_idx:

  Scalar row index.

- config:

  The SimulationConfig used to create the grid.

## Value

A named task specification list containing the task and grid indices,
data and fit specifications, task context, and RNG seed.
