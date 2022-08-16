library(bayesim)
library(extraDistr)
library(brms)

get_b <- function(mu, eta) {
  b <- (1 / mu) * log((-1 / eta) * log(1 / 2) + 1)
  return(b)
}
get_a <- function(mu, eta) {
  # a of extraDistr
  a <- get_b(mu, eta) * eta
  return(a)
}

n <- 10000   # number of testvalues

x = seq(from = 0 , to = 10 , length.out=n)

# PDF test
layout(matrix(1:3, ncol = 3))
plot(x, bayesim::dgompertz(x, mu = 2, eta = 0.1), type="l", ylab = "Density", main = "apex-after-origin Gompertz(mu=2, eta=0.1)")
plot(x, bayesim::dgompertz(x, mu = 1, eta = 1), type="l", ylab = "Density", main = "apex-at-origin Gompertz(mu=1, eta=1)")
plot(x, bayesim::dgompertz(x, mu = 2, eta = 3), type="l", ylab = "Density", main = "no-apex Gompertz(mu=2, eta=3)")

# PDF comparison to a reference
layout(matrix(1:3, ncol = 3))
plot(x, bayesim::dgompertz(x, mu = 2, eta = 0.1) - extraDistr::dgompertz(x, get_a(2, 0.1), get_b(2, 0.1)), ylab = "Density difference", main = "apex-after-origin Gompertz(mu=2, eta=0.1)")
plot(x, bayesim::dgompertz(x, mu = 1, eta = 1) - extraDistr::dgompertz(x, get_a(1, 1), get_b(1, 1)), ylab = "Density difference", main = "apex-at-origin Gompertz(mu=1, eta=1)")
plot(x, bayesim::dgompertz(x, mu = 2, eta = 3) - extraDistr::dgompertz(x, get_a(2, 3), get_b(2, 3)), ylab = "Density difference", main = "no-apex Gompertz(mu=2, eta=3)")

# Quantile test
x = seq(from = 0 , to = 1 , length.out=n)
layout(matrix(1:3, ncol = 3))
plot(x, bayesim::qgompertz(x, mu = 2, eta = 0.1), type="l", ylab = "Quantile", main = "apex-after-origin Gompertz(mu=2, eta=0.1)")
plot(x, bayesim::qgompertz(x, mu = 1, eta = 1), type="l", ylab = "Quantile", main = "apex-at-origin Gompertz(mu=1, eta=1)")
plot(x, bayesim::qgompertz(x, mu = 2, eta = 3), type="l", ylab = "Quantile", main = "no-apex Gompertz(mu=2, eta=3)")

# Quantile comparison to a reference
x = seq(from = 0 , to = 1 , length.out=n)
layout(matrix(1:3, ncol = 3))
plot(x, bayesim::qgompertz(x, mu = 2, eta = 0.1) - extraDistr::qgompertz(x, get_a(2, 0.1), get_b(2, 0.1)), ylab = "Quantile difference", main = "apex-after-origin Gompertz(mu=2, eta=0.1)")
plot(x, bayesim::qgompertz(x, mu = 1, eta = 1) - extraDistr::qgompertz(x, get_a(1, 1), get_b(1, 1)), ylab = "Quantile difference", main = "apex-at-origin Gompertz(mu=1, eta=1)")
plot(x, bayesim::qgompertz(x, mu = 2, eta = 3) - extraDistr::qgompertz(x, get_a(2, 3), get_b(2, 3)), ylab = "Quantile difference", main = "no-apex Gompertz(mu=2, eta=3)")


layout(matrix(1:3, ncol = 3))
y = bayesim::rgompertz(n, mu = 2, eta = 0.1)
hist(y, main = c(paste("Median:", median(y)), " for RNG of apex-after-origin Gompertz(mu=2, eta=0.1)"))
y = bayesim::rgompertz(n, mu = 1, eta = 1)
hist(y, main = c(paste("Median:", median(y)), " for RNG of apex-at-origin Gompertz(mu=1, eta=1)"))
y = bayesim::rgompertz(n, mu = 2, eta = 3)
hist(y, main = c(paste("Median:", median(y)), " for RNG of no-apex Gompertz(mu=2, eta=3)"))

# Now all the R code is tested, now test the BRMS family
n = 1000
a = rnorm(n)
data = list(a = a, y = bayesim::rgompertz(n, exp(0.5 * a + 1), 0.2))
layout(1)
hist(data$y)

fit1 <- brm(
  y ~ 1 + a,
  init = 0.1,
  data = data,
  family = bayesim::gompertz(),
  stanvars = bayesim::gompertz()$stanvars,
  backend = "cmdstanr",
  cores = 4
)

summary(fit1)

plot(fit1)

brms::pp_check(fit1)

brms::conditional_effects(fit1)

