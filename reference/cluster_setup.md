# Convenience Function to set up a cluster used for multiprocessing

Convenience Function to set up a cluster used for multiprocessing

## Usage

``` r
cluster_setup(ncores = 2, cluster_type = "FORK", debug = FALSE, outfile = NULL)
```

## Arguments

- ncores:

  Numbers of processes the cluster should use

- cluster_type:

  "FORK" or "PSOCK". Fork is faster but doesn't work on Windows

- debug:

  TRUE if the cluster log should be written to file

- outfile:

  Path where the cluster log is written to in debug mode.

## Value

A parallel cluster object registered with doParallel.
