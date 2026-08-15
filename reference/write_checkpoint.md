# Write checkpoint

Atomically writes a checkpoint containing the task grid (ledger) and
results. Uses a write-temp-then-rename protocol to ensure consistency.

## Usage

``` r
write_checkpoint(
  result_path,
  task_grid,
  task_results,
  config_fingerprint,
  checkpoint_format = "rds"
)
```

## Arguments

- result_path:

  Character string giving the base path for results. If NULL, the
  function returns NULL immediately (no checkpointing).

- task_grid:

  A tibble/data.frame containing the task grid with status information.
  Each row represents a task with columns for task_id, status, and other
  task metadata.

- task_results:

  A list of `bayesim_task_result` objects containing results from
  completed tasks.

- config_fingerprint:

  Character string containing a hash of the configuration for validation
  purposes.

- checkpoint_format:

  Character scalar naming the checkpoint storage format.

## Value

Invisible checkpoint ID (integer), or NULL if result_path is NULL.

## Details

The checkpoint directory structure is:


    checkpoints/
    +-- cp_000001/
        +-- meta.json         # checkpoint metadata
        +-- ledger.rds        # task grid with status
        +-- results.rds       # metrics + diagnostics per task
        +-- checksums.json    # file integrity checksums

Atomic write protocol:

1.  Create temporary directory with `.tmp` suffix

2.  Write all files to temporary directory

3.  Compute and write checksums

4.  Validate checksums in temporary directory

5.  Rename temporary directory to final name (atomic operation)

6.  Update latest.json pointer

If any step fails, the temporary directory is cleaned up and an error is
thrown. This ensures that partial writes never appear as valid
checkpoints.

## See also

[`read_checkpoint()`](https://sims1253.github.io/bayesim/reference/read_checkpoint.md),
[`init_checkpoint_dir()`](https://sims1253.github.io/bayesim/reference/init_checkpoint_dir.md)

## Examples

``` r
if (FALSE) { # \dontrun{
checkpoint_id <- write_checkpoint(
  result_path = "/path/to/results",
  task_grid = task_grid,
  task_results = task_results,
  config_fingerprint = "abc123hash"
)
} # }
```
