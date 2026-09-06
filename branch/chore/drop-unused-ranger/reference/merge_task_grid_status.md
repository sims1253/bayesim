# Merge Task Grid Status from Checkpoint

Updates a freshly-generated task grid with terminal status values from a
checkpoint. Pending tasks remain pending; only tasks that were completed
(success/failed/skipped) in the checkpoint are updated.

## Usage

``` r
merge_task_grid_status(fresh_grid, checkpoint_grid)
```

## Arguments

- fresh_grid:

  Task grid tibble from
  [`create_task_grid()`](https://sims1253.github.io/bayesim/reference/create_task_grid.md).

- checkpoint_grid:

  Task grid tibble from checkpoint.

## Value

Task grid tibble with status merged from checkpoint.

## Details

The function preserves the deterministic task ordering and RNG seeds
from the fresh grid, only updating status (and any recorded
`stop_reason`) for tasks that were terminal in the checkpoint. This
ensures resumed runs maintain identical reproducibility guarantees.
Policy-stopped rows return to pending with no stop reason; a resumed run
that stops before executing them re-marks them in
[`execute_tasks()`](https://sims1253.github.io/bayesim/reference/execute_tasks.md).
