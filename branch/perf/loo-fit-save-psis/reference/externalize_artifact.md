# Externalize large artifact

Moves a large artifact to an external RDS file and replaces the
in-memory object with a pointer containing metadata about the
externalized file.

## Usage

``` r
externalize_artifact(artifact, artifacts_dir, task_id, field_name)
```

## Arguments

- artifact:

  An R object to externalize

- artifacts_dir:

  Directory path where artifacts should be stored

- task_id:

  Task identifier used to construct the filename

- field_name:

  Name of the field being externalized (e.g., "draws", "fit")

## Value

A list with pointer information containing:

- external: TRUE (indicating this is an external pointer)

- path: Full path to the external file

- hash: Hash of the artifact for integrity checking

- size: Size of the artifact in bytes

## See also

[`write_rds_atomic()`](https://sims1253.github.io/bayesim/reference/write_rds_atomic.md),
[`compute_hash()`](https://sims1253.github.io/bayesim/reference/compute_hash.md)
