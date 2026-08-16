# Validate that resume does not promise artifacts that were never persisted.

Retention is runtime policy and therefore does not alter study identity,
but widening it after terminal outcomes exist cannot recreate discarded
draws, fits, predictions, data, diagnostics, or warnings. Reject that
ambiguity at the resume seam before any pending task is executed.

## Usage

``` r
validate_resume_retention(requested, persisted, checkpoint)
```
