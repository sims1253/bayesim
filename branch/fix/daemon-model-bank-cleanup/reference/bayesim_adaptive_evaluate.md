# Evaluate adaptive precision and return a persistable decision snapshot.

Evaluate adaptive precision and return a persistable decision snapshot.

## Usage

``` r
bayesim_adaptive_evaluate(
  task_results_so_far,
  task_grid,
  stop_on,
  config,
  completed_rounds = completed_replicate_rounds(task_grid),
  checked_at = sum(task_grid$status == "success", na.rm = TRUE)
)
```
