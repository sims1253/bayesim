# Forward-simulate the response for a single draw across fitters.

Dispatches by fitter class:

- LinearRegressionFitter (and fitters that expose fit_result\$draws):
  reconstruct a one-draw fit_result and call predict_fit(); since
  predict_fit consumes the ambient RNG state for noise, this is a single
  fresh Gaussian draw at the selected theta.

- BrmsFitter: call brms::posterior_predict(fit, newdata, draw_ids) for a
  single draw.

- Fallback: try predict_fit() with a sliced fit_result; error on
  failure.

## Usage

``` r
.forward_sim_y(fitter, fit_result, newdata, draw_id, theta_vec, voi, resp)
```

## Details

Returns a numeric vector of length nrow(newdata).
