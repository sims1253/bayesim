# Construct a fixed-truth data generator

Returns a generator closure that delegates to a user-supplied
`draw_data` function, pinning `data_bundle$true_params` to the supplied
`truth` for every task. Use this when the data-generating truth is known
and constant (e.g. a fixed effect size sweep).

## Usage

``` r
fixed_truth_generator(truth, draw_data)
```

## Arguments

- truth:

  Named numeric vector; the data-generating parameters, used as
  `true_params` and `vars_of_interest` for every generated bundle.

- draw_data:

  Function with signature `(data_spec, task_ctx)` returning a list with
  at least `train` (data.frame). May optionally return `response`,
  `test`, `meta`.

## Value

A generator function `(data_spec, task_ctx) -> data_bundle`.

## Details

The `draw_data` function must have signature
`(data_spec, task_ctx) -> list(train, test, response, ...)` WITHOUT
`true_params`/`vars_of_interest` — those are injected by the factory
from `truth`. It must consume the ambient RNG state (do not call
`set.seed`/[`withr::with_seed`](https://withr.r-lib.org/reference/with_seed.html)
inside).

## Examples

``` r
if (FALSE) { # \dontrun{
gen <- fixed_truth_generator(
  truth = c(beta = 1, sigma = 1),
  draw_data = function(data_spec, task_ctx) {
    n <- data_spec$n %||% 20L
    x <- stats::rnorm(n)
    y <- x + stats::rnorm(n)
    list(train = data.frame(y = y, x = x), response = "y")
  }
)
} # }
```
