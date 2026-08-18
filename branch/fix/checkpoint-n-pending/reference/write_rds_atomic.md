# Write RDS file atomically

Writes an RDS file using an atomic write operation to prevent corruption
from interrupted writes. Writes to a temporary file first, then renames.

## Usage

``` r
write_rds_atomic(x, path)
```

## Arguments

- x:

  Object to write to RDS

- path:

  Character string giving the file path

## Value

Invisible NULL. Called for side effect.
