# Convert SimulationConfig to Plain List for Hashing

Converts an S7 SimulationConfig object to a plain list suitable for
hashing or serialization. Excludes runtime-specific fields like
result_path and checkpoint_every.

## Usage

``` r
as_config_spec(config)
```

## Arguments

- config:

  An S7 SimulationConfig object.

## Value

A named list containing the configuration specification.

## Examples

``` r
if (FALSE) { # \dontrun{
config <- simulation_config(...)
spec <- as_config_spec(config)
# spec can now be serialized or hashed
} # }
```
