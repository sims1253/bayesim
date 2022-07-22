library(bayesim)
library(brms)

get_alpha <- function(mu, phi) {
  return(mu * phi)
}
get_beta <- function(mu, phi) {
  return((1 - mu) * phi)
}

n <- 10000   # number of testvalues

x = seq(from = 0.01, to = 0.99 , length.out=n)

# PDF test
layout(matrix(1:2, ncol = 2))
plot(x, bayesim::dbeta_custom(x, mu = 0.1, phi = 1), type="l", ylab = "Density", main = "left-leaning dbeta_custom(x, mu = 0.1, phi = 1)")
plot(x, bayesim::dbeta_custom(x, mu = 0.5, phi = 2), type="l", ylab = "Density", main = "uniform dbeta_custom(x, mu = 0.5, phi = 2)")


plot(x, bayesim::dbeta_custom(x, mu = 0.1, phi = 1) - dbeta(x, get_alpha(0.1, 1), get_beta(0.1, 1)), type="l", ylab = "Density", main = "left-leaning dbeta_custom(x, mu = 0.1, phi = 1)")
plot(x, bayesim::dbeta_custom(x, mu = 0.5, phi = 2) - dbeta(x, get_alpha(0.5, 2), get_beta(0.5, 2)), type="l", ylab = "Density", main = "uniform dbeta_custom(x, mu = 0.5, phi = 2)")

# RNG test
y = bayesim::rbeta_custom(n, mu = 0.1, phi = 1)
hist(y, main = c("RNG-plot ", paste("Mean:", mean(y)), " for RNG of left-leaning dbeta_custom(x, mu = 0.1, phi = 1)"))
y = bayesim::rbeta_custom(n, mu = 0.5, phi = 2)
hist(y, main = c("RNG-plot ", paste("Mean:", mean(y)), " for RNG of uniform dbeta_custom(x, mu = 0.5, phi = 2)"))

# Now all the R code is tested, now test the BRMS family
n = 1000
a = rnorm(n)
data = list(a = a, y = bayesim::rbeta_custom(n, bayesim::inv_logit(0.5 * a + 1), 2))
layout(1)
hist(data$y)

fit1 <- brm(
  y ~ 1 + a,
  data = data,
  family = bayesim:::beta_custom(),
  stanvars = bayesim:::beta_custom()$stanvars,
  backend = "cmdstan",
  cores = 4
)

summary(fit1)

plot(fit1)

brms::pp_check(fit1)

brms::conditional_effects(fit1)
