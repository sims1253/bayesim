# Write checksums for multiple files

Computes checksums for specified files in a directory and writes them to
a checksums.json file.

## Usage

``` r
write_checksums(dir_path, files)
```

## Arguments

- dir_path:

  Character string giving the directory path

- files:

  Character vector of file names (relative to dir_path)

## Value

Invisible NULL. Called for side effect.
