# Get latest valid checkpoint

Scans backwards from the latest checkpoint to find the most recent
checkpoint with valid checksums. Used for corruption recovery.

## Usage

``` r
get_latest_valid_checkpoint(result_path, config_fingerprint = NULL)
```

## Arguments

- result_path:

  Character string giving the base path for results.

- config_fingerprint:

  Optional character string. If provided, also validates that the
  checkpoint's fingerprint matches.

## Value

A checkpoint list (from
[`read_checkpoint()`](https://sims1253.github.io/bayesim/reference/read_checkpoint.md)),
or NULL if no valid checkpoint is found.

## Details

This function implements the corruption recovery algorithm:

1.  Get list of all checkpoint IDs

2.  Start from the highest ID

3.  Try to read each checkpoint (validates checksums)

4.  If fingerprint is provided, also validate fingerprint

5.  Return first valid checkpoint found

6.  Return NULL if no valid checkpoint found

## See also

[`read_checkpoint()`](https://sims1253.github.io/bayesim/reference/read_checkpoint.md),
[`list_checkpoints()`](https://sims1253.github.io/bayesim/reference/list_checkpoints.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Find latest valid checkpoint (any fingerprint)
checkpoint <- get_latest_valid_checkpoint("/path/to/results")

# Find latest valid checkpoint with matching fingerprint
checkpoint <- get_latest_valid_checkpoint(
  "/path/to/results",
  config_fingerprint = expected_fingerprint
)
} # }
```
