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
  [`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md)

- newdata:

  Data frame with new observations for prediction. If NULL, predictions
  are generated for the original training data.

- seed:

  Optional integer seed for reproducible predictions

## Value

A list containing:

- `predicted_mean`: Vector of mean predictions (length N)

- `predicted_samples`: Matrix of posterior predictive samples (S x N;
  draws as rows, observations as columns)

- `predicted_sd`: Vector of prediction standard deviations (length N)

- Additional fitter-specific outputs

`predicted_samples` follows the same orientation convention as
[`log_lik_matrix()`](https://sims1253.github.io/bayesim/reference/log_lik_matrix.md)
and
[`predict_epred()`](https://sims1253.github.io/bayesim/reference/predict_epred.md):
all matrices are draws x observations (S rows, N columns).
