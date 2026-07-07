# Check if Resumable Run Exists

Checks whether a valid checkpoint exists that can be resumed. A run can
be resumed if:

- Both run_manifest.json and latest.json exist

- latest.json points to a valid checkpoint_id

- The referenced checkpoint can be read and validated

## Usage

``` r
can_resume(result_path)
```

## Arguments

- result_path:

  Character; path to results directory containing checkpoints.

## Value

TRUE if a valid run can be resumed, FALSE otherwise.

## Examples

``` r
if (FALSE) { # \dontrun{
if (can_resume("/path/to/results")) {
  summary <- get_resume_summary("/path/to/results")
  cli::cli_alert_info("Found {summary$n_completed} completed tasks")
}
} # }
```
