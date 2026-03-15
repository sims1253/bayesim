# Lighten task result after checkpointing

Removes heavy objects from a task result that has been checkpointed.
This reduces memory usage while keeping lightweight summary data needed
for the final result construction.

## Usage

``` r
lighten_task_result(task_result, retain)
```

## Arguments

- task_result:

  A bayesim_task_result object

- retain:

  Character vector of retention options specifying what to keep

## Value

A lightened bayesim_task_result object with heavy fields removed. If all
heavy fields are removed and only minimal data remains, may return a
lightweight summary object instead.

## Details

After a task result has been checkpointed to disk, we can safely remove
heavy objects from memory since they can be reconstructed from the
checkpoint. This function:

- Always keeps: task_id, status, metrics, timing, error

- Keeps diagnostics only if "diagnostics" is in retain

- Removes: draws, predictions, fit, data (always, since checkpointed)

- Keeps warnings only if "warnings" is in retain

The rationale is that metrics and diagnostics are needed for the final
summary dataframe construction, while heavy objects (fit, draws, data)
can be loaded from checkpoint if needed for detailed analysis.

## See also

[`execute_tasks()`](https://sims1253.github.io/bayesim/reference/execute_tasks.md),
[`write_checkpoint()`](https://sims1253.github.io/bayesim/reference/write_checkpoint.md)
