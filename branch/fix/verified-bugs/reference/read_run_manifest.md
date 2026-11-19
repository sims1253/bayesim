# Read run manifest

Reads the run manifest from the result path.

## Usage

``` r
read_run_manifest(result_path)
```

## Arguments

- result_path:

  Character string giving the base path for results.

## Value

A list containing the run manifest, or NULL if not found.

## Details

The run manifest contains:

- `run_schema_version` - Schema version for the run structure

- `result_schema_version` - Schema version for result files

- `config_fingerprint` - Hash of the configuration

- `created` - Timestamp when the run was created
