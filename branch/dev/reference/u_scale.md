# Scale uniform samples to scaled ranks

Transforms uniform samples into scaled ranks for ECDF analysis. For a
matrix of samples, computes ranks within each column (chain).

## Usage

``` r
u_scale(x)
```

## Arguments

- x:

  A matrix of uniform samples with N rows and L columns (L chains)

## Value

A matrix of scaled ranks in the range 0 to 1
