# Initialize checkpoint directory

Creates the checkpoint directory structure for a new simulation run.
This includes the base result path, checkpoints subdirectory, run
manifest, and initial latest pointer.

## Usage

``` r
init_checkpoint_dir(
  result_path,
  config_fingerprint,
  config_spec = NULL,
  checkpoint_format = "rds"
)
```

## Arguments

- result_path:

  Character string giving the base path for results. If NULL, the
  function returns NULL immediately (no checkpointing).

- config_fingerprint:

  Character string containing a hash of the configuration. This is
  stored in the manifest for validation during resume operations.

- config_spec:

  Optional normalized configuration spec to persist in the run manifest
  for future rehydration.

- checkpoint_format:

  Character scalar naming the checkpoint storage format.

## Value

Invisible path to the result directory, or NULL if result_path is NULL.

## Details

The directory structure created is:


    result_path/
    +-- run_manifest.json    # run-level metadata and schema versions
    +-- latest.json          # pointer to latest valid checkpoint ID
    +-- checkpoints/         # directory for checkpoint snapshots

The run manifest contains:

- `run_schema_version` - Schema version for the run structure

- `result_schema_version` - Schema version for result files

- `config_fingerprint` - Hash of the configuration

- `created` - Timestamp when the run was created

## See also

[`write_checkpoint()`](https://sims1253.github.io/bayesim/reference/write_checkpoint.md),
[`read_checkpoint()`](https://sims1253.github.io/bayesim/reference/read_checkpoint.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result_path <- init_checkpoint_dir(
  result_path = "/path/to/results",
  config_fingerprint = "abc123hash"
)
} # }
```
