# Create a timer object

Creates a timer object with methods for tracking elapsed time across
execution phases.

## Usage

``` r
make_timer()
```

## Value

A list with the following methods:

- [`start()`](https://rdrr.io/r/stats/start.html) - Start or restart the
  timer

- [`stop()`](https://rdrr.io/r/base/stop.html) - Stop the timer

- `elapsed()` - Get elapsed time in seconds

## Examples

``` r
timer <- make_timer()
timer$start()
Sys.sleep(0.1)
timer$stop()
timer$elapsed()
#> [1] 0.1014268
```
