# Render a simulation-study report

Renders a standard Quarto HTML report for a `bayesim_simulation_result`,
covering the study design table,
[`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md),
parameter recovery plots per estimand, credible-interval coverage, and
SBC rank ECDF panels (when a rank metric was computed). The report
template lives at `inst/report/simulation-report.qmd`.

Each section is wrapped in `tryCatch`, so a missing metric (e.g. no rank
data, no recorded truths) degrades gracefully instead of failing the
whole render.

Requires the `quarto` R package AND the Quarto CLI. If the CLI is not
available, an informative error is thrown pointing to
<https://quarto.org>.

`render_report()` was previously named
[`report()`](https://sims1253.github.io/bayesim/reference/report.md);
the old name collided with the generic of the easystats *report* package
and now lives on as a deprecated alias (see
[`report()`](https://sims1253.github.io/bayesim/reference/report.md)).

## Usage

``` r
render_report(
  result,
  output_file = "bayesim-report.html",
  open = interactive()
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

## Value

The path to the rendered HTML file (invisibly).

## See also

[`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md),
[`plot_recovery()`](https://sims1253.github.io/bayesim/reference/plot_recovery.md),
[`plot_coverage()`](https://sims1253.github.io/bayesim/reference/plot_coverage.md),
[`plot_rank_ecdf()`](https://sims1253.github.io/bayesim/reference/plot_rank_ecdf.md).

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_simulation(config, progress = FALSE)
render_report(result, output_file = "my-study.html")
} # }
```
