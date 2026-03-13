# Validate schema compatibility

Checks whether the schema versions in a run manifest are compatible with
the current package version.

## Usage

``` r
validate_schema_compatibility(manifest)
```

## Arguments

- manifest:

  A run manifest list from
  [`read_run_manifest()`](https://sims1253.github.io/bayesim/reference/read_run_manifest.md).

## Value

Logical. TRUE if schemas are compatible, FALSE otherwise.

## Details

Currently, this function checks for exact version matches. In the
future, it may support backward compatibility for certain version
ranges.
