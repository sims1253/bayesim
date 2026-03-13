# Custom LOO object builder

Builds a loo object that contains any pointwise criterion, acting as
elpd for compatibility.

## Usage

``` r
custom_loo_object(pointwise_criterion, psis_object = NULL)
```

## Arguments

- pointwise_criterion:

  vector of criterion values for each observation

- psis_object:

  PSIS object for psis diagnostics

## Value

a loo object, containing a criterion, disguised as elpd

## Details

LOO currently has hardcoded expectations of elpd as part of loo objects
so to use loo objects, we have to disguise other criterions as elpd.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a custom loo object from pointwise criterion values
custom_loo_object(c(-1.2, -0.8, -1.5))
} # }
```
