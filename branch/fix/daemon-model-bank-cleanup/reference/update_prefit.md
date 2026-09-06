# Update a prefit with new data (model bank path)

Runs [`stats::update()`](https://rdrr.io/r/stats/update.html) with
`recompile = FALSE`. BEFORE updating, verifies that the task data +
formula + family produce the same Stan *data structure* as the prefit
(via
[`brms::make_standata()`](https://paulbuerkner.com/brms/reference/standata.html)).
brms does NOT warn when `recompile = FALSE` is passed against
structurally incompatible data — it silently reuses the compiled binary
against a wrong model frame — so this explicit structural comparison is
the only reliable guard. Raises a fatal
[`bayesim_internal_error()`](https://sims1253.github.io/bayesim/reference/bayesim_internal_error.md)
on any mismatch (e.g. different predictor count, new factor levels,
missing variables).

## Usage

``` r
update_prefit(
  prefit,
  fitter,
  data_bundle,
  seed,
  formula,
  family,
  prior,
  stanvars
)
```

## Details

Note: we compare the Stan *data structure* (`make_standata` field names
and the design-matrix column count `K`), NOT `make_stancode`, because
the latter embeds data-derived default-prior values (e.g. the
intercept's `student_t` location = `mean(y)`) that vary across datasets
without affecting whether the compiled binary is reusable.
