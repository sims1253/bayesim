# Run Schema Version

Version identifier for checkpoint format compatibility. Increment this
when the on-disk checkpoint format changes in a way that breaks backward
compatibility.

Increments on any incompatible on-disk format change to the run-level
structure (directory layout, manifest format, etc.).

## Usage

``` r
RUN_SCHEMA_VERSION

RUN_SCHEMA_VERSION
```

## Format

An object of class `integer` of length 1.

An object of class `integer` of length 1.
