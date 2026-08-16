# Prune old checkpoint commit directories

Prune old checkpoint commit directories

## Usage

``` r
prune_checkpoints(result_path, keep = 2L)
```

## Arguments

- result_path:

  Checkpoint result directory.

- keep:

  Number of newest checkpoint commits to keep; `Inf` keeps all. Removes
  commit directories only; immutable outcome shards and ledger history
  are never pruned.

## Value

Invisible vector of removed checkpoint IDs.
