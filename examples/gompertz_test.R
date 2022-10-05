library(bayesim)
library(extraDistr)
library(brms)

# get_b <- function(mu, b) {
#   b <- (1 / mu) * log((-1 / b) * log(1 / 2) + 1)
#   return(b)
# }
# get_a <- function(mu, b) {
#   # a of extraDistr
#   a <- get_b(mu, b) * b
#   return(a)
# }

get_a <- function(mu, beta) {
  # a of extraDistr
  a <- -(beta * log(0.5)) / (exp(mu * beta) - 1)
  return(a)
}

n <- 10000   # number of testvalues

x = seq(from = 0.01, to = 10 , length.out=n)

# PDF test
layout(matrix(1:3, ncol = 3))
plot(x, bayesim::dgompertz(x, mu = 2, b = 0.1), type="l", ylab = "Density", main = "Gompertz(mu=2, b=0.1)")
plot(x, bayesim::dgompertz(x, mu = 1, b = 1), type="l", ylab = "Density", main = "Gompertz(mu=1, b=1)")
plot(x, bayesim::dgompertz(x, mu = 2, b = 3), type="l", ylab = "Density", main = "Gompertz(mu=2, b=3)")

# PDF comparison to a reference
layout(matrix(1:3, ncol = 3))
plot(x, bayesim::dgompertz(x, mu = 2, b = 0.1) - extraDistr::dgompertz(x, get_a(2, 0.1), 0.1), ylab = "Density difference", main = "Gompertz(mu=2, b=0.1)")
plot(x, bayesim::dgompertz(x, mu = 1, b = 1) - extraDistr::dgompertz(x, get_a(1, 1), 1), ylab = "Density difference", main = "Gompertz(mu=1, b=1)")
plot(x, bayesim::dgompertz(x, mu = 2, b = 3) - extraDistr::dgompertz(x, get_a(2, 3), 3), ylab = "Density difference", main = "Gompertz(mu=2, b=3)")

# Quantile test
x = seq(from = 0.01, to = 0.99, length.out=n)
layout(matrix(1:3, ncol = 3))
plot(x, bayesim::qgompertz(x, mu = 2, b = 0.1), type="l", ylab = "Quantile", main = "Gompertz(mu=2, b=0.1)")
plot(x, bayesim::qgompertz(x, mu = 1, b = 1), type="l", ylab = "Quantile", main = "Gompertz(mu=1, b=1)")
plot(x, bayesim::qgompertz(x, mu = 2, b = 3), type="l", ylab = "Quantile", main = "Gompertz(mu=2, b=3)")

# Quantile comparison to a reference
x = seq(from = 0.01, to = 0.99, length.out=n)
layout(matrix(1:3, ncol = 3))
plot(x, bayesim::qgompertz(x, mu = 2, b = 0.1) - extraDistr::qgompertz(x, get_a(2, 0.1), 0.1), ylab = "Quantile difference", main = "Gompertz(mu=2, b=0.1)")
plot(x, bayesim::qgompertz(x, mu = 1, b = 1) - extraDistr::qgompertz(x, get_a(1, 1), 1), ylab = "Quantile difference", main = "Gompertz(mu=1, b=1)")
plot(x, bayesim::qgompertz(x, mu = 2, b = 3) - extraDistr::qgompertz(x, get_a(2, 3), 3), ylab = "Quantile difference", main = "Gompertz(mu=2, b=3)")


layout(matrix(1:3, ncol = 3))
y = bayesim::rgompertz(n, mu = 2, b = 0.1)
hist(y, main = c(paste("Median:", median(y)), " for RNG Gompertz(mu=2, b=0.1)"))
y = bayesim::rgompertz(n, mu = 1, b = 1)
hist(y, main = c(paste("Median:", median(y)), " for RNG Gompertz(mu=1, b=1)"))
y = bayesim::rgompertz(n, mu = 2, b = 3)
hist(y, main = c(paste("Median:", median(y)), " for RNG Gompertz(mu=2, b=3)"))

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

