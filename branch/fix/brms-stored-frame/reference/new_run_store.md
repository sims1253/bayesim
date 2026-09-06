# Create a run store backed by memory or the filesystem.

Create a run store backed by memory or the filesystem.

## Usage

``` r
new_run_store(
  result_path = NULL,
  config_fingerprint = NULL,
  config_spec = NULL,
  checkpoint_format = "rds",
  keep_checkpoints = 2L,
  retention_spec = NULL,
  run_policy_spec = NULL
)
```

## Arguments

- result_path:

  NULL for memory, otherwise a filesystem directory.

- config_fingerprint:

  Study fingerprint.

- config_spec:

  Optional manifest specification.

- checkpoint_format:

  Checkpoint serialization format.

- keep_checkpoints:

  Number of checkpoint commit directories to retain. Pruning removes old
  commit directories only; immutable outcome shards and ledger history
  are never pruned, so durable storage grows roughly linearly with
  completed tasks.

## Value

An internal run-store object with initialize/read/write methods.
