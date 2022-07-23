library(bayesim)
library(VGAM)
library(testthat)

n <- 100000 # number of testvalues
eps <- 1e-3
x <- seq(from = eps, to = 1 - eps, length.out = n)

n_small <- 10
mus <- seq(from = eps, to = 1 - eps, length.out = n_small)
ps <- seq(from = eps, to = 10, length.out = n_small)

# mus_r <- seq(from = 1 + eps, to = 10, length.out = n_small)
ps_r <- seq(from = 0.6 + eps, to = 10, length.out = n_small)
accepted_means_eps <- 0.02
p_acceptable_failures <- 0.05

get_q <- function(mu, p) {
  value <- -log(2) / log1p(-mu^p)
  return(value)
}

# check many shape parameters on pdf and qdf
for (p in ps) {
  for (m in mus) {
    expect_eps(bayesim::dkumaraswamy(x, mu = m, p = p), VGAM::dkumar(x, p, get_q(m, p)), eps)
    #expect_eps(bayesim::qkumaraswamy(x, mu = m, p = p), VGAM::qkumar(x, p, get_q(m, p)), eps)
    diff <- abs(bayesim::qkumaraswamy(x, mu = m, p = p) - VGAM::qkumar(x, p, get_q(m, p)))
    if(isTRUE(any(diff > eps))) {
      print(c("quantiles different: p, m", p, m))
      data <- VGAM::qkumar(x, p, get_q(m, p))
      layout(matrix(1:3, ncol = 3))
      plot(x, bayesim::qkumaraswamy(x, mu = m, p = p))
      plot(x, VGAM::qkumar(x, p, get_q(m, p)))
      plot(x, bayesim::qkumaraswamy(x, mu = m, p = p) - VGAM::qkumar(x, p, get_q(m, p)))
      invisible(readline(prompt="Press [enter] to continue"))
    }
  }
}

