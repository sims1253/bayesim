# Clean old checkpoints

Removes old checkpoints, keeping only the N most recent ones.

## Usage

``` r
clean_old_checkpoints(result_path, keep_n = 5)
```

## Arguments

- result_path:

  Character string giving the base path for results.

- keep_n:

  Integer. Number of most recent checkpoints to keep. Default is 5.

## Value

Invisible number of checkpoints removed.

## Details

This function can be used to manage disk space by removing old
checkpoints that are no longer needed. It always keeps the latest
checkpoint and any checkpoints that are more recent than the keep_n
threshold.
