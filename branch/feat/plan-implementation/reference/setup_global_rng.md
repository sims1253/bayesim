# Set up global RNG for simulation

Sets RNG kind to L'Ecuyer-CMRG and initializes with global seed. Must be
called once at simulation start.

## Usage

``` r
setup_global_rng(seed)
```

## Arguments

- seed:

  Integer seed for the simulation

## Value

The initial `.Random.seed` state (invisibly)

## Examples

``` r
if (FALSE) { # \dontrun{
setup_global_rng(42)
} # }
```
