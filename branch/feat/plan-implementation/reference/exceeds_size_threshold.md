# Check if result row exceeds size threshold

Determines whether a task result exceeds a specified memory threshold,
useful for deciding whether to externalize large artifacts.

## Usage

``` r
exceeds_size_threshold(task_result, threshold_bytes = 5 * 1024 * 1024)
```

## Arguments

- task_result:

  A bayesim_task_result object

- threshold_bytes:

  Maximum allowed size in bytes. Default is 5 MB.

## Value

TRUE if the task_result exceeds the threshold, FALSE otherwise

## Examples

``` r
if (FALSE) { # \dontrun{
result <- list(metrics = list(rmse = 0.1), diagnostics = list(rhat = 1.01))
exceeds_size_threshold(result)
exceeds_size_threshold(result, threshold_bytes = 10)
} # }
```
