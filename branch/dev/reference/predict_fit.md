# Generate Predictions

Generate predictions from a fitted model.

## Usage

``` r
predict_fit(fitter, fit_result, newdata = NULL, seed = NULL)
```

## Arguments

- fitter:

  An S7 Fitter object

- fit_result:

  A `bayesim_fit_result` object from
  [`fit()`](https://sims1253.github.io/bayesim/reference/fit.md)

- newdata:

  Data frame with new observations for prediction. If NULL, predictions
  are generated for the original training data.

- seed:

  Optional integer seed for reproducible predictions

## Value

A list containing:

- `predicted_mean`: Vector of mean predictions (N)

- `predicted_samples`: Matrix of posterior predictive samples (N x S)

- `predicted_sd`: Vector of prediction standard deviations (N)

- Additional fitter-specific outputs
