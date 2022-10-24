library(bayesim)

eps <- 1e-12 # 2 digits more, than sim
# intervals
unit_int <- c(eps, 1 - eps)
mu_unit_int <- c(0.1, 0.9)
pos_int <- c(eps, 200)
shape_int <- c(0.1, 20)
n_x <- 1000
n_p <- 20

density_lookup_generator(
  n_param = n_p,
  n_x = n_x,
  mu_int = mu_unit_int,
  mu_link = cauchit,
  shape_int = shape_int,
  x_int_prelink = unit_int,
  density_fun = dcauchitnormal,
  density_name = "cauchitnormal"
)
density_lookup_generator(
  n_param = n_p,
  n_x = n_x,
  mu_int = mu_unit_int,
  mu_link = cloglog,
  shape_int = shape_int,
  x_int_prelink = unit_int,
  density_fun = dcloglognormal,
  density_name = "cloglognormal"
)
density_lookup_generator(
  n_param = n_p,
  n_x = n_x,
  mu_int = mu_unit_int,
  mu_link = logit,
  shape_int = shape_int,
  x_int_prelink = unit_int,
  density_fun = dlogitnormal,
  density_name = "logitnormal"
)
density_lookup_generator(
  n_param = n_p,
  n_x = n_x,
  mu_int = pos_int,
  shape_int = shape_int,
  x_int_prelink = pos_int,
  x_link = exp,
  density_fun = dlognormal_custom,
  density_name = "lognormal"
)
density_lookup_generator(
  n_param = n_p,
  n_x = n_x,
  mu_int = pos_int,
  shape_int = shape_int,
  x_int_prelink = pos_int,
  x_link = exp,
  density_fun = dsoftplusnormal, density_name = "softplusnormal"
)
