# Prune old checkpoint snapshots

Prune old checkpoint snapshots

## Usage

``` r
prune_checkpoints(result_path, keep = 2L)
```

## Arguments

- result_path:

  Checkpoint result directory.

- keep:

  Number of newest complete snapshots to keep; `Inf` keeps all.

## Value

Invisible vector of removed checkpoint IDs.
