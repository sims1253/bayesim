# List available checkpoints

Returns a vector of all checkpoint IDs available in the result path.

## Usage

``` r
list_checkpoints(result_path)
```

## Arguments

- result_path:

  Character string giving the base path for results.

## Value

Integer vector of checkpoint IDs, sorted in ascending order. Returns an
empty integer vector if no checkpoints exist or the checkpoints
directory doesn't exist.

## Details

This function scans the checkpoints subdirectory and extracts IDs from
directory names matching the pattern `cp_XXXXXX`.

## See also

[`read_checkpoint()`](https://sims1253.github.io/bayesim/reference/read_checkpoint.md),
[`get_latest_valid_checkpoint()`](https://sims1253.github.io/bayesim/reference/get_latest_valid_checkpoint.md)

## Examples

``` r
if (FALSE) { # \dontrun{
checkpoint_ids <- list_checkpoints("/path/to/results")
# Returns: c(1L, 2L, 3L, 5L, 10L)
} # }
```
