// Bayesian linear regression for the CmdStanFitter worked example (D2).
// Mirrors LinearRegressionFitter's model: y ~ normal(X * beta, sigma).
// log_lik (pointwise) and mu (epred) generated quantities for elpd / prediction
// metrics. log_lik and mu are vector[N]; posterior reads them as S x N.
data {
  int<lower=1> N;            // observations
  int<lower=1> K;            // predictors (excluding intercept)
  matrix[N, K + 1] X;        // design matrix, first column is the intercept
  vector[N] y;
}
parameters {
  vector[K + 1] beta;        // coefficients (intercept first)
  real<lower=0> sigma;
}
model {
  beta ~ normal(0, 10);
  sigma ~ normal(0, 10);
  y ~ normal(X * beta, sigma);
}
generated quantities {
  vector[N] log_lik;
  vector[N] mu;              // posterior expectation of the mean (epred)
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | X[n] * beta, sigma);
    mu[n] = X[n] * beta;
  }
}
