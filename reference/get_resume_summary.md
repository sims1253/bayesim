# Get Resume Summary

Returns a summary of what would be resumed from a checkpoint. Useful for
informing users about resume state before actually loading and resuming.

## Usage

``` r
get_resume_summary(result_path)
```

## Arguments

- result_path:

  Character; path to results directory containing checkpoints.

## Value

A list with elements:

- `checkpoint_id`: ID of the checkpoint

- `n_total`: Total number of tasks

- `n_completed`: Number of completed (success + failed) tasks

- `n_pending`: Number of pending tasks

Returns NULL if no valid resume state exists.

## Examples

``` r
if (FALSE) { # \dontrun{
summary <- get_resume_summary("/path/to/results")
if (!is.null(summary)) {
  cli::cli_alert_info("Checkpoint: {summary$checkpoint_id}")
  cli::cli_alert_info("Completed: {summary$n_completed}/{summary$n_total}")
  cli::cli_alert_info("Remaining: {summary$n_pending}")
}
} # }
```
