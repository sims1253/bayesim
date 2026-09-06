# Build brms control list from stan_args

Maps the user-facing `stan_args` fields (`adapt_delta`, `max_treedepth`)
into the `control` list that brms/Stan accept, alongside the `init` and
`threads` passthroughs. Returns a named list of extra arguments suitable
for splicing into
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) /
[`stats::update()`](https://rdrr.io/r/stats/update.html) via
[`do.call()`](https://rdrr.io/r/base/do.call.html), or an empty list
when `stan_args` is NULL.

## Usage

``` r
build_stan_args_call(stan_args)
```

## Arguments

- stan_args:

  Named list or NULL.

## Value

A named list (possibly empty) of brms arguments.
