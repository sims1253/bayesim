# Resolve a family specification to a canonical brms family object

Accepts a brms family object, a base `family` object, a string (resolved
via
[`brms::brmsfamily()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)),
or NULL (defaults to `brms::gaussian()`). Base `family` objects (e.g.
from [`stats::gaussian()`](https://rdrr.io/r/stats/family.html)) are
canonicalized through
[`brms::brmsfamily()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
so that [`gaussian()`](https://rdrr.io/r/stats/family.html) and
`brmsfamily("gaussian")` produce identical objects (and thus identical
model-spec hashes and Stan code).

## Usage

``` r
resolve_family(family)
```

## Arguments

- family:

  A brms/base family object, a string, or NULL.

## Value

A brms family object.
