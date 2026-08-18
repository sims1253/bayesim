# Verify checksums for files in a directory

Reads the checksums.json file and verifies all files match their
recorded checksums.

## Usage

``` r
verify_checksums(dir_path)
```

## Arguments

- dir_path:

  Character string giving the directory path

## Value

Logical. TRUE if all checksums match, FALSE otherwise.
