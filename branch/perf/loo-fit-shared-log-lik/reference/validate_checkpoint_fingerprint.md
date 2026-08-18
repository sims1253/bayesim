# Validate checkpoint fingerprint

Checks whether a checkpoint's configuration fingerprint matches the
expected fingerprint. This is used during resume to ensure the
checkpoint is compatible with the current configuration.

## Usage

``` r
validate_checkpoint_fingerprint(checkpoint, config_fingerprint)
```

## Arguments

- checkpoint:

  A checkpoint list returned by
  [`read_checkpoint()`](https://sims1253.github.io/bayesim/reference/read_checkpoint.md).

- config_fingerprint:

  Character string containing the expected configuration fingerprint.

## Value

Logical. TRUE if fingerprints match exactly, FALSE otherwise.

## Details

The config fingerprint is a hash of the normalized configuration that
excludes runtime-only fields (result_path, checkpoint_every, progress).
A mismatch indicates that the checkpoint was created with a different
configuration and should not be used for resuming.

## See also

[`read_checkpoint()`](https://sims1253.github.io/bayesim/reference/read_checkpoint.md)

## Examples

``` r
if (FALSE) { # \dontrun{
checkpoint <- read_checkpoint("/path/to/results")
if (validate_checkpoint_fingerprint(checkpoint, expected_fingerprint)) {
  # Safe to resume from this checkpoint
} else {
  # Configuration has changed, cannot resume
}
} # }
```
