# Deprecated alias for render_report()

`report()` was renamed to
[`render_report()`](https://sims1253.github.io/bayesim/reference/render_report.md)
because the old name collided with the generic of the easystats *report*
package. The alias forwards to
[`render_report()`](https://sims1253.github.io/bayesim/reference/render_report.md)
and emits a deprecation warning once per session; it will be removed in
a future release.

## Usage

``` r
report(
  result,
  output_file = "bayesim-report.html",
  open = interactive(),
  estimands = NULL
)
```

## Arguments

- result:

  A `bayesim_simulation_result` from
  [`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md).

- output_file:

  Path to the rendered HTML output file (default
  `"bayesim-report.html"`).

- open:

  Logical scalar. When `TRUE` (the default in interactive sessions) the
  rendered report is opened in a viewer/browser. Forwarded to
  [`quarto::quarto_render()`](https://quarto-dev.github.io/quarto-r/reference/quarto_render.html)
  via its `open` handling where supported; on systems without that
  argument the file path is still returned.

- estimands:

  Ignored. The argument was accepted (and documented as informational
  only) by `report()` but never used; it is kept in the signature purely
  so existing calls do not error.
