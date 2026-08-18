# Coerce a formula to a brmsformula

Plain formulas are wrapped via `brms::::brmsformula()` equivalent
([`brms::bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html));
existing brmsformulas are returned unchanged.

## Usage

``` r
resolve_formula(formula)
```

## Arguments

- formula:

  A formula or brmsformula.

## Value

A brmsformula object.
