library(bayesim)
library(extraDistr)
library(brms)

get_lambda <- function(mu, alpha) {
  lambda <- mu * (alpha - 1)
  # lambda of extraDistr is inverted, compared to lambda in my pdf
  return(1/lambda)
}

n <- 10000   # number of testvalues

x = seq(from = 0, to = 10 , length.out=n)

# PDF test
layout(matrix(1:2, ncol = 2))
plot(x, bayesim::dlomax(x, mu = 2, alpha = 2), type="l", ylab = "Density", main = "high starting Lomax(mu=2, alpha=2)")
plot(x, bayesim::dlomax(x, mu = 8, alpha = 2), type="l", ylab = "Density", main = "low starting Lomax(mu=8, alpha=2)")

# PDF extraDistr test (extraDistr::dlomax(x = 0) == 0?)
#layout(matrix(1:2, ncol = 2))
#plot(x, extraDistr::dlomax(x, get_lambda(2, 2), 2), type="l", ylab = "extraDistr::Density", main = "high starting Lomax(mu=2, alpha=2)")
#plot(x, extraDistr::dlomax(x, get_lambda(8, 2), 2), type="l", ylab = "extraDistr::Density", main = "low starting Lomax(mu=8, alpha=2)")


x = seq(from = 10^-6, to = 10 , length.out=n)
#PDF comparison
layout(matrix(1:2, ncol = 2))
plot(x, bayesim::dlomax(x, mu = 2, alpha = 2) - extraDistr::dlomax(x, get_lambda(2, 2), 2), ylab = "Density difference", main = "high starting Lomax(mu=2, alpha=2)")
plot(x, bayesim::dlomax(x, mu = 8, alpha = 2) - extraDistr::dlomax(x, get_lambda(8, 2), 2), ylab = "Density difference", main = "low starting Lomax(mu=8, alpha=2)")

x = seq(from = 0 , to = 1 , length.out=n)
# Quantile test
layout(matrix(1:2, ncol = 2))
plot(x, bayesim::qlomax(x, mu = 2, alpha = 2), type="l", ylab = "Quantile", main = "high starting Lomax(mu=2, alpha=2)")
plot(x, bayesim::qlomax(x, mu = 8, alpha = 2), type="l", ylab = "Quantile", main = "low starting Lomax(mu=8, alpha=2)")

# Quantile difference
layout(matrix(1:2, ncol = 2))
plot(x, bayesim::qlomax(x, mu = 2, alpha = 2) - extraDistr::qlomax(x, get_lambda(2, 2), 2), ylab = "Quantile difference", main = "high starting Lomax(mu=2, alpha=2)")
plot(x, bayesim::qlomax(x, mu = 8, alpha = 2) - extraDistr::qlomax(x, get_lambda(8, 2), 2), ylab = "Quantile difference", main = "low starting Lomax(mu=8, alpha=2)")

# RNG test
y = bayesim::rlomax(n, mu = 2, alpha = 2)
hist(log(y), main = c("log-plot ", paste("Mean:", mean(y)), " for RNG of high starting Lomax(mu=2, alpha=2)"))
y = bayesim::rlomax(n, mu = 8, alpha = 2)
hist(log(y), main = c("log-plot ", paste("Mean:", mean(y)), " for RNG of low starting Lomax(mu=8, alpha=2)"))

# Now all the R code is tested, now test the BRMS family
n = 1000
a = rnorm(n)
data = list(a = a, y = bayesim::rlomax(n, exp(0.5 * a + 1), 2))
layout(1)
hist(log(data$y))

fit1 <- brm(
  y ~ 1 + a,
  data = data,
  family = bayesim::lomax(),
  stanvars = bayesim::lomax()$stanvars,
  backend = "cmdstan",
  cores = 4
)

summary(fit1)

plot(fit1)

brms::pp_check(fit1)

brms::conditional_effects(fit1)
