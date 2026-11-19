# Get next checkpoint ID

Determines the next available checkpoint ID by examining existing
checkpoint directories.

## Usage

``` r
get_next_checkpoint_id(result_path)
```

## Arguments

- result_path:

  Character string giving the base path for results.

## Value

Integer. The next checkpoint ID (1 if no checkpoints exist).

## Details

Checkpoint directories are named with the format `cp_XXXXXX` where X
represents a zero-padded six-digit number. This function scans existing
directories and returns the maximum ID + 1.
