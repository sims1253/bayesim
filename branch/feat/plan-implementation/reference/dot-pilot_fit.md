# Run the one-time preconditioning pilot fit via fit_model().

Uses the fitter's fit_model() with a fixed seed so the preconditioning
fit is reproducible across tasks (the worker's ambient RNG stream is
independent of the seed passed here, which only governs the MCMC/NIG
draw generation for the pilot). Returns the bayesim_fit_result.

## Usage

``` r
.pilot_fit(fitter, pilot_bundle, fit_spec)
```
