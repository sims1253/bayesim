library(BBmisc)       # to surpress unnecassary BRMS/Stan output
library(posterior)    # used to extract BRMS data
library(brms)         # to test custom BRMS families

# used only for an internal function

#' Check that |a - b| < eps. Works with scalars and vectors on any input.
#' For vector input, one may defined, how many |a - b| >= eps are acceptable with r argument.
#'
#' @param a numeric scalar or vector a to be compared
#' @param b numeric scalar or vector b to be compared
#' @param eps numeric scalar or vector, setting the max differences, eps > 0
#' @param r optional numeric scalar (r = 0 in defualt), relative number of values, that may have an difference > eps.
#' Used internally, tlo calculate acceptable amount of deviances. Calculated absolute value will be floored.
#'
#' @md
#' @details For vector/scalar combinations, allowed are (r has to always be a scalar!):
#'   * all scalar
#'   * a vector, b scalar, eps scalar -> compare a values to be close enough to b
#'   * a scalar, b vector, eps scalar -> compare b values to be close enough to a
#'   * a scalar, b scalar, eps vector -> compare the difference a, b too all eps
#'   -> case not intended, but would work
#'   * a vector, b scalar, eps vector -> compare each vector-difference entry against each eps
#'   * all vectors -> each entry of |a - b| is compared to the same entry in eps
#'   -> different vector lengths != 1 dissallowed!expect_brms_family
#'
#' @return success or failure with message
#'
#' @export
#'
#' @examples library(testthat)
#' library(bayesim)
#' expect_eps(1, 1.1, 0.2) # should pass
#' expect_eps(-1, -1.1, 0.2) # should pass
#' expect_eps(-1.00001, -1.00001, 1e-4) # should pass
#' expect_eps(c(0, 1, 2), c(1, 1, 2), 1e-4, 0.4) # should pass (only 1/3 were wrong, but 40% was allowed)
#' expect_eps(c(0, 1, 3), c(1, 1, 2), 1e-4, 1/3) # should fail (2/3 were wrong, but only 1/3 was allowed)
#' expect_eps(1, 1.3, 0.2) # should fail (|a - b| = 0.3 > 0.2)
#' expect_eps(c(1, 1.1), c(1, 1.1, 1.2), 0.5) # should produce an error (different vector-lengths)
#' expect_eps(1, 1.1, -0.2) # should produce an error (negative eps)
#' expect_eps(c(0, 1, 2), c(1, 1, 2), 1e-4, -0.4) # should produce an error (r too small)
#' expect_eps(c(0, 1, 2), c(1, 1, 2), 1e-4, 1.4) # should produce an error (r too big)
expect_eps <- function(a, b, eps, r = 0.0) {
  # first check, that the type is always a scalar or vector
  all_numeric <- is.numeric(a) && is.numeric(b) && is.numeric(eps) && is.numeric(r)
  if (isFALSE(all_numeric)) {
    stop("At least one of the vectors was not a scalar/vector.")
  }

  # then check, that r is only a scalar
  if (isTRUE(length(r) > 1)) {
    stop("The relative number of tolerated deviances r has to be a scalar.")
  }

  # also check, that r is in range [0, 1]
  if (isTRUE(r < 0 || r > 1)) {
    stop("The relative number of tolerated deviances r has to be in [0, 1].")
  }

  # then check, that all eps are >= 0
  # (changed to >= given the r might also prove interesting for i.e. integer comparisons)
  if (isTRUE(any(eps < 0))) {
    stop("Tried checking against negative differences.")
  }

  # then check, if vectors are of same length
  vector_length <- max(length(a), length(b), length(eps))
  tolerated_deviances <- floor(vector_length * r)

  a_wrong <- length(a) > 1 && length(a) != vector_length
  b_wrong <- length(b) > 1 && length(b) != vector_length
  eps_wrong <- length(eps) > 1 && length(eps) != vector_length
  if (a_wrong || b_wrong || eps_wrong) {
    stop("Used different length of vectors in test.")
  }

  # last check, the actual value difference
  eps_comparison_wrong <- abs(a - b) > eps
  number_deviances <- sum(eps_comparison_wrong, na.rm = TRUE)
  if (isTRUE(number_deviances <= tolerated_deviances)) {
    succeed()
  } else {
    if (isTRUE(tolerated_deviances > 0)) {
      message <- paste("In expect_eps: ", toString(number_deviances), " of ", toString(vector_length), " were bigger, then eps. Only ", toString(tolerated_deviances), " allowed!", sep = "")
    } else {
      message <- paste("In expect_eps: ", toString(number_deviances), " of ", toString(vector_length), " were bigger, then eps. None were allowed!", sep = "")
    }
    fail(message)
  }
}


#' Distribution RNG test vingette
#'
#' @param rng_fun RNG function under test
#' @param metric_mu Metric to be used on RNG data (usually mean or median)
#' @param mus Metric data used as RNG argument and to be compared to (usually mean or median)
#' @param shapes Shape arguemnts
#' @param mu_eps Acceptable difference of |mu - metric_mu(rng_fun)
#' @param p_acceptable_failures Acceptable rate of failure, relative value of difference bigger mu_eps
#' @param mu_link Default=identity, optional link-function argument, for example
#' useful in link-normal-distributions
#'
#' @return Nothing actually, just wraps the test
#' @export
#'
#' @examples library(testthat)
#' library(bayesim)
#' mus <- seq(from = 1 + eps, to = 20, length.out = 10)
#' phis <- seq(from = 2 + eps, to = 20, length.out = 10)
#' test_rng(rng_fun=bayesim::rbetaprime, metric_mu=mean, n=10000, mus=mus,
#'          shapes=phis, mu_eps=0.2, p_acceptable_failures=0.05)
test_rng <- function(rng_fun, metric_mu, n, mus, shapes, mu_eps, p_acceptable_failures, mu_link=identity) {
  if(isFALSE(is.function(rng_fun) && is.function(metric_mu) && is.function(mu_link))) {
    stop("RNG-, Metric- or mu_link-function argument was not a function!")
  }
  if(isFALSE(is.numeric(n) && length(n) == 1 && n >= 1)) {
    stop("n must be a numeric, positive scalar!")
  }
  if(isFALSE(length(mu_eps) == 1)) {
    stop("mu_eps has to be a scalar in this test-function!")
  }

  lenx <- length(shapes)
  leny <- length(mus)
  expected_mus <- rep(mus, times=leny)
  rng_mus <- vector(length=lenx * leny)

  for (x in 1:lenx) {
    for (y in 1:leny) {
      rng_mus[(x-1) * leny + y] <- mu_link(metric_mu(rng_fun(n, mu = mus[y], shapes[x])))
    }
  }
  expect_eps(rng_mus, expected_mus, mu_eps, p_acceptable_failures)
}

#' Check, if a is bigger than b. expect_smaller(a, b) can be done,
#' by expect_bigger(b, a) with swapped arguments.
#'
#' For vector/scalar combinations, allowed are:
#' both scalar
#' a scalar, b vector -> compare a to all b values
#' a vector, b scalar -> compare all a values to the threshold b
#' both vectors -> compare each entry to the other entry
#' -> different lengths != 1 dissallowed!
#'
#' @param a number to be tested, should be real scalar or vector
#' @param b tested against, should be real scalar or vector
#'
#' @return expect(a > threshold)
#'
#' @export
#'
#' @examples expect_bigger(1, 0) # should pass
#' expect_bigger(0, 1) # should fail
#' expect_bigger(c(1, 2), 0) # should pass
#' expect_bigger(c(0, 2), 1) # should fail
#' expect_bigger(c(1, 2), c(0, 3)) # should fail
#' expect_bigger(c(1, 2, 3), c(0, 1)) # should produce an error
expect_bigger <- function(a, b) {

  # check, that both are numerical data
  if (!is.numeric(a) || !is.numeric(b)) {
    stop("All inputs have to be numeric data.")
  }

  # check, if vectors both are vectors, both are of equal length
  vector_length <- max(length(a), length(b))
  a_wrong <- length(a) > 1 && length(a) != vector_length
  b_wrong <- length(b) > 1 && length(b) != vector_length
  if (a_wrong || b_wrong) {
    stop("Used different lengths of vectors in test")
  }

  # now compare all a and b values
  all_bigger <- all(a > b)
  if (isTRUE(all_bigger)) {
    succeed()
  } else {
    fail("At least one number a was smaller than b")
  }
}

# TODO: does this need doc+test??? (Don't think so, given it is just a wrapper)
expect_brms_family <- function(n_data_sampels=1000, ba=0.5, intercept=1, shape=2, link, family, rng, shape_name, thresh=0.05, num_parallel=2, debug=FALSE, data_threshold=NULL, seed=1337) {
  if(isFALSE(is.character(shape_name) && length(shape_name) == 1)) {
    stop("The shape name argument has to be a single string")
  }
  posterior_fit <- construct_brms(n_data_sampels, ba, intercept, shape, link, family, rng, num_parallel, seed=seed, suppress_output=!debug, data_threshold=data_threshold)
  #plot(posterior_fit)

  ba_succ  <- test_brms_quantile(posterior_fit, "b_a", ba, thresh, debug)
  int_succ <- test_brms_quantile(posterior_fit, "b_Intercept", intercept, thresh, debug)
  sha_succ <- test_brms_quantile(posterior_fit, shape_name, shape, thresh, debug)
  success  <- ba_succ && int_succ && sha_succ

  if(debug && !success) {
    print("Data were not recovered correctly! Print plot.")
    plot(posterior_fit)
  }

  if(success)
    succeed()
  else
    fail("One or more variables have not been recovered correctly! You may set debug true and check the plot.")
}

# TODO: does this need doc+test??? (Don't think so, given it is just a wrapper)
construct_brms <- function(n_data_sampels, ba, intercept, shape, link, family, rng, num_parallel, seed=NULL, suppress_output=TRUE, data_threshold=NULL) {
  if(isFALSE(is.function(family) && is.function(rng))) {
    stop("family or rng argument were not a function!")
  }

  old_seed <- .Random.seed
  set.seed(seed)
  a <- rnorm(n_data_sampels)
  y_data <- rng(n_data_sampels, link(ba * a + intercept), shape)
  if(!is.null(data_threshold))
    y_data <- limit_data(y_data, data_threshold)
  set.seed(old_seed)

  data <- list(a = a, y = y_data)

  if(isTRUE(suppress_output)) {
    # if printout suppression is wished, use suppressAll as wrapper
    BBmisc::suppressAll({
      posterior_fit <- brms::brm(
        y ~ 1 + a,
        data = data,
        family = family(),
        stanvars = family()$stanvars,
        backend = "cmdstanr",
        chains = num_parallel,
        cores = num_parallel,
        #seed = seed, #no seed for BRMS. BRMS RNG code is to be tested to, so would be not sensible
        silent = 2,
        refresh = 0,
        init = 0.1
      )
    })
  }
  else {
    # and if not, do nothing special
    posterior_fit <- brms::brm(
      y ~ 1 + a,
      data = data,
      family = family(),
      stanvars = family()$stanvars,
      backend = "cmdstanr",
      chains = num_parallel,
      cores = num_parallel,
      #seed = seed, #no seed for BRMS. BRMS RNG code is to be tested to, so would be not sensible
      silent = 2,
      refresh = 0,
      init = 0.1
    )
  }


  return(posterior_fit)
}

limit_data <- function(data, limits) {
  # check that the limit is usable
  if(length(limits) != 2 && limits[1] > limits[2]) {
    stop("If the limits is to be used, it has to be of size 2.
         Index 1 is the lower bound and 2 is the upper bound. Specify NA, if a bound is unused.")
  }

  # if so, use the applicable limit (If one uses to na, well. What are you trying to achieve? :)
  if(!is.na(limits[1]))
    data[data < limits[1]] <- limits[1]
  if(!is.na(limits[2]))
    data[data > limits[2]] <- limits[2]

  # now return the data
  return(data)

}

#' Check, that data of the posterior is close enough to the reference data.
#'
#' @param posterior_data Data fitted and drawn
#' @param name Name of the variable to check
#' @param reference Reference value to check against
#' @param thresh Scalar or 2-length vector of quantile bounds.
#' For scalar constructs bound as [tresh, 1-thresh]
#'
#' @return succeess, fail or error
#' @export
#'
#' @examples n <- 1000
#' a <- rnorm(n)
#' ba_in <- 0.5
#' data <- list(a = a, y = bayesim::rbetaprime(n, exp(ba_in * a + 1), 2))
#' fit1 <- brm(y ~ 1 + a, data = data, family = bayesim::betaprime(),
#'   stanvars = bayesim::betaprime()$stanvars, backend = "cmdstan",
#'   cores = 4, silent = 2, refresh = 0)
#' test_brms_quantile(fit1, "b_a", ba_in, 0.025)
#' TODO: how to test this function (and is this to be tested)?
test_brms_quantile <- function(posterior_data, name, reference, thresh, debug = FALSE) {
  if(isTRUE(length(thresh) == 1))
    bounds = c(thresh, 1-thresh)
  else if(isTRUE(length(thresh) == 2))
    bounds = thresh
  else
    stop("The quantile-thresholds can only be a scalar, or a 2 length vector")

  calculated <- posterior::extract_variable_matrix(posterior_data, variable = name)
  quantiles <- unname(quantile(calculated, probs = bounds))
  if(debug) {
    print(paste("At: ", bounds, " the quantile is: ", quantiles))
    print(paste("The supplied reference value was: ", reference))
  }

  # if(quantiles[1] < reference && reference < quantiles[2])
  #   succeed()
  # else
  #   fail("The reference value was not within the quantiles of the given posterior data!")
  return(quantiles[1] < reference && reference < quantiles[2])
}


#' Vector combinator. Used in a few link-normal tests (to simplify indexing, hehe).
#' Not meant to be used by users. (Hence no at-export)
#'
#' @param a Vector a
#' @param b Vector b
#'
#' @return combination c(c(a[1], b[1]), c(a[1], b[2]), ...)
#'
#' @examples print(bayesim:::vector_combinator(seq(0, 1, length.out=3), seq(0, 2, length.out=3)))
vector_combinator <- function(a, b) {
  # TODO: not only is this unused by now. Also I do think, this might be possible in vanilla R.
  length_a <- length(a)
  length_b <- length(b)
  total_length <- length_a * length_b
  output <- list(length = total_length)
  for (x in 1:length_a) {
    for (y in 1:length_b) {
      output[[(x-1) * length_a + y]] <- c(a[x], b[y])
    }
  }
  return(output)
}

#' Yannick's super dooper density lookup data generator!
#' This is meant for development only. Use at your own risk (and sanity). You have been warned!
#'
#' @param n_param Number of paramerter for each mu and shape. Will yield n_param^2 permutations of both.
#' @param n_x Number of the x values.
#' @param eps Allowed precision tolerance.
#' @param mu_int Interval of the densities mu variable.
#' @param shape_int Interval of the densities shape variable.
#' @param density_fun The PDF to be saved.
#' @param density_name The name of the whole pre-calculated thingy. Probably the one of the density.
#' @param x_int_prelink Interval of x-values, before the x_link function is used on them.
#' @param x_link Link function to be used on x_int_prelink. Default is identity.
#'
#' @return Nothing, but saves input and reference files for later use!
#' @export
#'
#' @examples library(bayesim)
#' eps <- 1e-6
#' density_lookup_generator  (mu_int=c(eps, 1-eps), shape_int=c(2,10), x_int_prelink=c(eps, 1-eps),
#'                           density_fun=bayesim::dcauchitnormal, density_name="cauchitnormal_demodata")
density_lookup_generator <- function(n_param=10, n_x=100, eps=10^-6, mu_int, shape_int, x_int_prelink, x_link=identity, density_fun, density_name) {
  # assert parameters (why all the assertions, if only I ever get to use it?)
  # dunno, just felt like it :P
  if(isFALSE(is.function(density_fun))) {
    stop("The given density_fun was no function. Expect density-function as input.")
  }
  if(isFALSE(isSingleString(density_name))) {
    stop("The given density_name has to be a single string.")
  }
  if(isFALSE(is.function(x_link))) {
    stop("The given x_link was no function. Expect link-function as input.")
  }
  if(isFALSE(isNum_len(mu_int, 2))) {
    stop("The mu_int interval has to be a two long vector")
  }
  if(isFALSE(isNum_len(shape_int, 2))) {
    stop("The shape_int interval has to be a two long vector")
  }
  if(isFALSE(isNum_len(x_int_prelink, 2))) {
    stop("The x_int_prelink interval has to be a two long vector")
  }
  if(isFALSE(isInt_len(n_param) && n_param > 0)) {
    stop("The n_param has to be a single whole integer bigger 0")
  }
  if(isFALSE(isInt_len(n_x) && n_x > 0)) {
    stop("The n_x has to be a single whole integer bigger 0")
  }
  if(isFALSE(isNum_len(eps) && eps > 0)) {
    stop("The eps has to be a scalar bigger 0")
  }

  # generate input data
  mus <- seq(from=mu_int[1], to=mu_int[2], length.out=n_param)
  shapes <- seq(from=shape_int[1], to=shape_int[2], length.out=n_param)

  x_ref <- x_link(seq(from=x_int_prelink[1], to=x_int_prelink[2], length.out=n_x))
  y_ref <- matrix(data=list(), nrow=n_param, ncol=n_param)

  # generate reference data
  for(outer in 1:n_param) {
    for(inner in 1:n_param) {
      mu <- mus[outer]
      shape <- shapes[inner]
      y_ref[[outer, inner]] <- density_fun(x_ref, mu, shape)
    }
  }

  inp_scalars <- data.frame(n=n_x, n_small=n_param, eps=eps)
  inp_vectors <- data.frame(mus=mus, shapes=shapes)
  x_data <- data.frame(x=x_ref)
  y_data <- data.frame(y=y_ref)

  return_data <- c(inp_scalars, inp_vectors, x_data)
  saveRDS(return_data, paste("tests/testthat/precalc_values/", density_name, "_refdata", sep=""))
  saveRDS(y_data, paste("tests/testthat/precalc_values/", density_name, "_refpdf", sep=""))
  # after reading the help, the issue with saveRDS (and save for that matter) is,
  # that it might be used incorrectly, which will not occur, if any bayesim dev implements them

  print(paste("Generated lookup data \"", density_name, "\" and saved to file in precalc_values folder for tests", sep=""))
}

isNum_len <- function(num, len=1) {
  return(is.numeric(num) && length(num) == len)
}

isInt_len <- function(int, len=1) {
  return((int %% 1 == 0) && length(int) == len)
}

isSingleString <- function(input) {
  return (is.character(input) && length(input) == 1)
}
