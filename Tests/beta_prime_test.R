library(bayesim)
library(extraDistr)
library(brms)

get_alpha <- function(mu, phi) {
  alpha <- mu * (phi + 1)
}

get_beta <- function(phi){
  beta <- phi + 2
}

n <- 10000   # number of testvalues

x = seq(from = 0.00001 , to = 10 , length.out=n)

layout(matrix(1:3, ncol = 3))
plot(x, bayesim::dbetaprime(x, mu = 1, phi = 2), type="l", ylab = "Density", main = "left-leaning Beta-Prime(mu=1,phi=2)")
plot(x, bayesim::dbetaprime(x, mu = 1, phi = 4), type="l", ylab = "Density", main = "rising-quick Beta-Prime(mu=1,phi=4)")
plot(x, bayesim::dbetaprime(x, mu = 8, phi = 2), type="l", ylab = "Density", main = "rising-slow Beta-Prime(mu=8,phi=2)")


layout(matrix(1:3, ncol = 3))
plot(x, bayesim::dbetaprime(x, mu = 1, phi = 2) - extraDistr::dbetapr(x, get_alpha(1, 2), get_beta(2)), ylab = "Difference", main = "left-leaning Beta-Prime(mu=1,beta=2)")
plot(x, bayesim::dbetaprime(x, mu = 1, phi = 4) - extraDistr::dbetapr(x, get_alpha(1, 4), get_beta(4)), ylab = "Difference", main = "rising-quick Beta-Prime(mu=1,beta=4)")
plot(x, bayesim::dbetaprime(x, mu = 8, phi = 2) - extraDistr::dbetapr(x, get_alpha(8, 2), get_beta(2)), ylab = "Difference", main = "rising-slow Beta-Prime(mu=8,beta=2)")



x = seq(from = 0 , to = 1 , length.out=n)
layout(matrix(1:3, ncol = 3))
plot(x, bayesim::qbetaprime(x, mu = 1, phi = 2), type="l", ylab = "Quantile", main = "left-leaning Beta-Prime(mu=1,phi=2)")
plot(x, bayesim::qbetaprime(x, mu = 1, phi = 4), type="l", ylab = "Quantile", main = "rising-quick Beta-Prime(mu=1,phi=4)")
plot(x, bayesim::qbetaprime(x, mu = 8, phi = 2), type="l", ylab = "Quantile", main = "rising-slow Beta-Prime(mu=8,phi=2)")

layout(matrix(1:3, ncol = 3))
plot(x, bayesim::qbetaprime(x, mu = 1, phi = 2) - extraDistr:::qbetapr(x, get_alpha(1, 2), get_beta(2)), ylab = "Quantile difference", main = "left-leaning Beta-Prime(mu=1,phi=2)")
plot(x, bayesim::qbetaprime(x, mu = 1, phi = 4) - extraDistr:::qbetapr(x, get_alpha(1, 4), get_beta(4)), ylab = "Quantile difference", main = "rising-quick Beta-Prime(mu=1,phi=4)")
plot(x, bayesim::qbetaprime(x, mu = 8, phi = 2) - extraDistr:::qbetapr(x, get_alpha(8, 2), get_beta(2)), ylab = "Quantile difference", main = "rising-slow Beta-Prime(mu=8,phi=2)")

layout(matrix(1:3, ncol = 3))
y = bayesim::rbetaprime(n, mu = 1, phi = 2)
hist(log(y), main = c(paste("Mean:", mean(y)), " for log(RNG) of left-leaning Beta-Prime(mu=1,phi=2)"))
y = bayesim::rbetaprime(n, mu = 1, phi = 4)
hist(log(y), main = c(paste("Mean:", mean(y)), " for log(RNG) of rising-quick Beta-Prime(mu=1,phi=4)"))
y = bayesim::rbetaprime(n, mu = 8, phi = 2)
hist(log(y), main = c(paste("Mean:", mean(y)), " for log(RNG) of rising-slow Beta-Prime(mu=8,phi=2)"))

# Now all the R code is tested, now test the BRMS family
n = 1000
a = rnorm(n)
data = list(a = a, y = bayesim::rbetaprime(n, exp(0.5 * a + 1), 2))
layout(1)
hist(data$y)

fit1 <- brm(
  y ~ 1 + a,
  data = data,
  family = bayesim::betaprime(),
  stanvars = bayesim::betaprime()$stanvars,
  backend = "cmdstan",
  cores = 4
)

summary(fit1)
# Errors seem a bit high, but nothing real suspicous

plot(fit1)
#look like beta-prime distributions. Should intercept and a look comparable?

brms::pp_check(fit1)
#may look fine, though the number scales are not liked by the plot

brms::conditional_effects(fit1)
