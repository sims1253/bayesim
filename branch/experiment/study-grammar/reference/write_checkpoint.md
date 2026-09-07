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
  checkpoint_format = "rds",
  keep_checkpoints = Inf,
  prior_results_df = NULL,
  prior_task_results = NULL,
  adaptive_next_check = NULL,
  adaptive_state = NULL,
  run_policy_spec = NULL,
  prior_checkpoint = NULL,
  return_state = FALSE,
  delta_store = FALSE
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

- keep_checkpoints:

  Positive integer; number of checkpoint commit directories to retain.
  Older commit directories are pruned only after the new checkpoint
  commit and `latest.json` have been written successfully. Pruning never
  removes immutable outcome shards or ledger history, so durable storage
  grows roughly linearly with completed tasks. `Inf` disables pruning
  for direct/internal callers.

- prior_results_df:

  Optional cached data frame of results from before the current
  execution. When supplied, it replaces the legacy read of the previous
  checkpoint.

- prior_task_results:

  Optional canonical task outcomes from before the current execution,
  used when migrating legacy checkpoints.

- adaptive_next_check:

  Optional persisted adaptive-check threshold.

- adaptive_state:

  Optional persistable adaptive precision snapshot.

- run_policy_spec:

  Optional serialized effective run policy.

- prior_checkpoint:

  Optional validated in-memory checkpoint state. The filesystem store
  supplies this to keep live append work proportional to newly completed
  tasks.

- return_state:

  Logical; return the validated checkpoint state instead of only its
  numeric ID.

- delta_store:

  Logical; write immutable outcome and ledger deltas instead of legacy
  full snapshots.

## Value

Invisible checkpoint ID (integer), or NULL if result_path is NULL.

## Details

The checkpoint directory structure is:


    checkpoints/
    +-- cp_000001/
        +-- meta.json         # checkpoint metadata
        +-- ledger.rds        # ledger view (delta in delta-store mode)
        +-- results.rds       # results view (delta in delta-store mode)
        +-- checksums.json    # file integrity checksums

With `delta_store = TRUE`, each checkpoint commit directory is an atomic
commit record: its `meta.json` plus the delta views above. Newly
completed outcomes are appended once as immutable, redundantly mirrored
shards under `outcomes/`, and task statuses under `ledger/` are a base
ledger plus status deltas. Readers reconstruct the accumulated run state
from the shards a commit references; they never rewrite history. With
`delta_store = FALSE` (legacy snapshot mode for direct/internal
callers), the checkpoint directory itself carries the snapshot views.

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
