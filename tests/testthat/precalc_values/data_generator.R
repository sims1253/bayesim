library(bayesim)

eps <- 1e-6
unit_int <- c(eps, 1-eps)
pos_int <- c(eps, 200)
shape_int <- c(2,10)

bayesim::density_lookup_generator(mu_int=unit_int, shape_int=shape_int, x_int_prelink=unit_int,
                                  density_fun=bayesim::dcauchitnormal, density_name="cauchitnormal")
bayesim::density_lookup_generator(mu_int=unit_int, shape_int=shape_int, x_int_prelink=unit_int,
                                  density_fun=bayesim::dcloglognormal, density_name="cloglognormal")
bayesim::density_lookup_generator(mu_int=unit_int, shape_int=shape_int, x_int_prelink=unit_int,
                                  density_fun=bayesim::dlogitnormal, density_name="logitnormal")
bayesim::density_lookup_generator(mu_int=unit_int, shape_int=shape_int, x_int_prelink=pos_int, x_link=exp,
                                  density_fun=bayesim::dlognormal_custom, density_name="lognormal")
bayesim::density_lookup_generator(mu_int=unit_int, shape_int=shape_int, x_int_prelink=pos_int, x_link=exp,
                                  density_fun=bayesim::dsoftplusnormal, density_name="softplusnormal")
