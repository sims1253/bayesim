# Compute the linear predictor X theta for a design data.frame and a draw.

coef_names names the additive coefficient columns; the intercept is
named "Intercept" and contributes a constant 1 (since predictor designs
do not include it as a column). Remaining coef_names must match columns
of newdata. Returns a length-nrow(newdata) numeric vector.

## Usage

``` r
.linear_predictor(newdata, theta_vec, coef_names, resp)
```
