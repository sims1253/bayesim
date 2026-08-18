# Write JSON file atomically

Writes a JSON file using an atomic write operation to prevent corruption
from interrupted writes. Writes to a temporary file first, then renames.

## Usage

``` r
write_json_atomic(x, path)
```

## Arguments

- x:

  Object to write to JSON

- path:

  Character string giving the file path

## Value

Invisible NULL. Called for side effect.

## Details

The function writes to a `.tmp` suffix file first, then uses
[`file.rename()`](https://rdrr.io/r/base/files.html) to atomically move
it to the target path. If the rename fails, a `bayesim_checkpoint_error`
is thrown.
