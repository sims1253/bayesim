
# used only for an internal function

#' Check that |a - b| < eps. Works with scalars and vectors on any input.
#' For vector input, one may defined, how many |a - b| >= eps are acceptable with r argument.
#' Given it is a "workhorse" used in almost all test-files, this function is also very well tested.
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
#' expect_eps(c(0, 1, 3), c(1, 1, 2), 1e-4, 1 / 3) # should fail (2/3 were wrong, but only 1/3 was allowed)
#' expect_eps(1, 1.3, 0.2) # should fail (|a - b| = 0.3 > 0.2)
#' expect_eps(c(1, 1.1), c(1, 1.1, 1.2), 0.5) # should produce an error (different vector-lengths)
#' expect_eps(1, 1.1, -0.2) # should produce an error (negative eps)
#' expect_eps(c(0, 1, 2), c(1, 1, 2), 1e-4, -0.4) # should produce an error (r too small)
#' expect_eps(c(0, 1, 2), c(1, 1, 2), 1e-4, 1.4) # should produce an error (r too big)
#' expect_eps(NA, 1, 0.1) # should produce an error, NAs are dissallowed
expect_eps <- function(a, b, eps, r = 0.0) {

  # then check, that r is only a scalar. Also check, that r is in range [0, 1)
  if (isFALSE(isNum_len(r) && r >= 0 && r < 1)) {
    stop("The relative number of tolerated deviances r has to be a scalar in [0, 1).")
  }

  # then check, that all eps are >= 0
  # (changed to >= given the r might also prove interesting for i.e. integer comparisons)
  if (isTRUE(any(eps < 0))) {
    stop("Tried checking against negative differences.")
  }

  # then check, if vectors are of same length
  vector_length <- max(length(a), length(b), length(eps))
  a_correct <- isNum_len(a) || isNum_len(a, vector_length)
  b_correct <- isNum_len(b) || isNum_len(b, vector_length)
  eps_correct <- isNum_len(eps) || isNum_len(eps, vector_length)
  # Now this construction should check is.numeric and should be immune to NA, I hope.
  if (!(a_correct && b_correct && eps_correct)) {
    stop("Used different length of numeric vectors in test. (Or vectors containing NAs)")
  }

  # For the test, calculate the absolute value, of how many entries in |a - b|
  # may be > eps. This is especially important for RNG tests, which will not always be precise
  # (obviously). Opposed to this in PDF comparisons, this figure is usually 0, as they
  # should always produce comparable results within a few eps machine precision.
  tolerated_deviances <- floor(vector_length * r)

  # last check, the actual value difference
  eps_comparison_wrong <- (abs(a - b) > eps)
  # convert the logical vector in a sum of how many entries were wrong
  number_deviances <- sum(eps_comparison_wrong, na.rm = TRUE)

  # at the end, check only the allowed number of entries were wrong
  if (isTRUE(number_deviances <= tolerated_deviances)) {
    succeed()
  } else {
    # in case of a failure, print how many entries have been incorrect and how many were allowed.
    if (isTRUE(tolerated_deviances > 0)) {
      message <- paste0("In expect_eps: ", toString(number_deviances), " of ", toString(vector_length), " were bigger, then eps. Only ", toString(tolerated_deviances), " allowed!")
    } else {
      message <- paste0("In expect_eps: ", toString(number_deviances), " of ", toString(vector_length), " were bigger, then eps. None were allowed!")
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
#' test_rng(
#'   rng_fun = bayesim::rbetaprime, metric_mu = mean, n = 10000, mus = mus,
#'   shapes = phis, mu_eps = 0.2, p_acceptable_failures = 0.05
#' )
test_rng <- function(rng_fun, metric_mu, n, mus, shapes, mu_eps, p_acceptable_failures, mu_link = identity) {

  # check, that all function arguments are actually functions
  # TODO: (Is it possible, to check, if they also take the correct arguments?)
  if (isFALSE(is.function(rng_fun) && is.function(metric_mu) && is.function(mu_link))) {
    stop("RNG-, Metric- or mu_link-function argument was not a function!")
  }

  # check the number of samples to be generated and cecked
  if (!(isInt_len(n) && n >= 1)) {
    stop("n must be an integer, positive scalar!")
  }

  # check, that compare eps is a scalar
  # all used moments should only deviate by eps in most tests)
  if (isFALSE(length(mu_eps) == 1)) {
    stop("mu_eps has to be a scalar in this test-function!")
  }

  # prepare the data, use a vector for ease of use
  # allows re-using the expect_eps.
  # As opposed to using a matrix, which would just complicate implementation and comparison.
  lenx <- length(shapes)
  leny <- length(mus)
  expected_mus <- rep(mus, times = leny)
  rng_mus <- vector(length = lenx * leny)

  # calculate rng data
  for (x in 1:lenx) {
    for (y in 1:leny) {
      rng_mus[(x - 1) * leny + y] <- mu_link(metric_mu(rng_fun(n, mu = mus[y], shapes[x])))
    }
  }
  # now the data was written, compare it
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

  # check, if vectors both are vectors, both are of equal length
  vector_length <- max(length(a), length(b))
  a_correct <- isNum_len(a) || isNum_len(a, vector_length)
  b_correct <- isNum_len(b) || isNum_len(b, vector_length)
  if (!(a_correct && b_correct)) {
    stop("Used different lengths of vectors in test, or non-numeric data")
  }

  # now compare all a and b values
  all_bigger <- all(a > b)
  if (isTRUE(all_bigger)) {
    succeed()
  } else {
    fail("At least one number a was smaller than b")
  }
}

#' BRMS family expect recovery. Tries linear baysian model y ~ 1 + a.
#' Checking of the arguments done in construct_brms.
#'
#' @param n_data_sampels How many samples per chain. Positive integer scalar. Default = 1000.
#' @param ba Linear data argument, real scalar. Default = 0.5.
#' @param intercept Intercept data argument, real scalar. Default = 1.
#' @param shape Shape argument of each distribution. Default = 2.
#' @param link Link function pointer used data. For positive bounded uses exp as example. No default.
#' @param family BRMS family under test.
#' @param rng function pointer of bespoke RNG for the family to be tested.
#' @param shape_name BRMS string of shape argument name. Single string.
#' @param num_parallel Chains calculated in parallel. Positive scalar integer. Default = 2.
#' @param seed Seed argument, so that input data is always the same in each test.
#' BRMS test does not test RNG and is not guarateed to fit on all data. Positive Integer scalar, Default = 1337.
#' Seed is stored before test and restored after it finished. If wants not to use a seed set to NA.
#' @param data_threshold Usually unused. But in rare cases, data too close at the boundary may cause trouble.
#' If so, set a two entry real vector c(lower, uppper). If one of them is NA, the data will not be capped for that boundary.
#' Default = Null, will be in R terms "invisible" and will not cap any input data.
#' @param thresh Acceptable threshold for quantiles of recovered arguments.
#' Scalar or 2-entry real vector within (0, 1).
#' Vector is used as is, scalar will be interpreted as c(tresh, 1-thresh).
#' Default = 0.05
#' @param debug Scalar Boolean argument, whether debug info is printed or not. Default = False.
#' Note: Will supress everything, besides errors (so tests stay clean).
#'
#' @return None
#' @export
#'
#' @examples library(testthat)
#' result <- expect_brms_family(link = exp, family = bayesim::betaprime, rng = bayesim::rbetaprime, shape_name = "phi")
#' print(result)
expect_brms_family <- function(n_data_sampels = 1000, ba = 0.5, intercept = 1, shape = 2, link, family, rng, shape_name, num_parallel = 2, seed = 1337, data_threshold = NULL, thresh = 0.05, debug = FALSE) {
  if (!isSingleString(shape_name)) {
    stop("The shape name argument has to be a single string")
  }
  posterior_fit <- construct_brms(n_data_sampels, ba, intercept, shape, link, family, rng, num_parallel, seed = seed, suppress_output = !debug, data_threshold = data_threshold)
  # plot(posterior_fit)

  ba_succ <- test_brms_quantile(posterior_fit, "b_a", ba, thresh, debug)
  int_succ <- test_brms_quantile(posterior_fit, "b_Intercept", intercept, thresh, debug)
  sha_succ <- test_brms_quantile(posterior_fit, shape_name, shape, thresh, debug)
  success <- ba_succ && int_succ && sha_succ

  if (debug && !success) {
    print("Data were not recovered correctly! Print plot.")
    fam_name <- posterior_fit$family$name
    debug <- paste0(
      fam_name, " expect_brms_family failed with inputs b_a=", ba,
      ", intercept=", intercept, " and shape=", shape
    )
    plot(posterior_fit, main = debug)
  }

  if (success) {
    succeed()
  } else {
    fail("One or more variables have not been recovered correctly! You may set debug true and check the plot.")
  }
}

#' Construct BRMS family for simple linear y ~ 1 + a model.
#'
#' @param n_data_sampels How many samples per chain. Positive integer scalar.
#' @param ba Linear data argument, real scalar.
#' @param intercept Intercept data argument, real scalar.
#' @param shape Shape argument of each distribution.
#' @param link Link function pointer used data. For positive bounded uses exp as example.
#' @param family BRMS family under test.
#' @param rng function pointer of bespoke RNG for the family to be tested.
#' @param shape_name BRMS string of shape argument name. Single string.
#' @param num_parallel Chains calculated in parallel. Positive scalar integer.
#' @param seed Seed argument, so that input data is always the same in each test.
#' BRMS test does not test RNG and is not guarateed to fit on all data.
#' Positive Integer scalar, Default = NA will do nothing. Seed is stored before and restored after.
#' @param data_threshold Usually unused. But in rare cases, data too close at the boundary may cause trouble.
#' If so, set a two entry real vector c(lower, uppper). If one of them is NA, the data will not be capped for that boundary.
#' Default = Null, will be in R terms "invisible" and will not cap any input data.
#' @param suppress_output Scalar Boolean argument. Default = TRUE supresses all prints.
#' Only exceptions will be printed. For testing reasons, to not spam the test-window.
#'
#' @return BRMS model for the specified family.
#' @export
#'
#' @examples posterior_fit <- construct_brms(1000, 0.5, 1.0, 2.0, exp, bayesim::betaprime, bayesim::rbetaprime, 2)
#' plot(posterior_fit)
#' # Only used internally with wrapper expect_brms_family, but I can't stop you anyways!
construct_brms <- function(n_data_sampels, ba, intercept, shape, link, family, rng,
                           num_parallel, seed = NA, data_threshold = NULL, suppress_output = TRUE) {
  if (!(is.function(family) && is.function(rng) && is.function(link))) {
    stop("family, rng or link argument were not a function!")
  }
  if (!(isInt_len(n_data_sampels) && n_data_sampels > 0)) {
    stop("n_data_sampels has to be a positive integer scalar")
  }
  if (!isNum_len(ba)) {
    stop("ba argument has to be a real scalar")
  }
  if (!isNum_len(intercept)) {
    stop("intercept argument has to be a real scalar")
  }
  if (!isNum_len(shape)) {
    stop("shape argument has to be a real scalar")
  }
  if (!(isNum_len(num_parallel) && num_parallel > 0)) {
    stop("num_prallel argument has to be a positive scalar integer")
  }
  if (!(isNum_len(seed) || is.na(seed))) {
    stop("seed argument if used has to be a real scalar. Else it is let default as NULL,
         which will not change the current RNG seed")
  }
  if (!isLogic_len(suppress_output)) {
    stop("The argument suppress_output has to be a single boolean")
  }


  if(!is.na(seed))
  {
    old_seed <- .Random.seed
    set.seed(seed)
  }

  a <- rnorm(n_data_sampels)
  y_data <- rng(n_data_sampels, link(ba * a + intercept), shape)
  if (!is.null(data_threshold)) {
    y_data <- limit_data(y_data, data_threshold)
  }

  if(!is.na(seed))
    set.seed(old_seed)


  data <- list(a = a, y = y_data)

  if (isTRUE(suppress_output)) {
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
        # seed = seed, #no seed for BRMS. BRMS RNG code is to be tested to, so would be not sensible
        silent = 2,
        refresh = 0,
        init = 0.1
      )
    })
  } else {
    # and if not, do nothing special
    posterior_fit <- brms::brm(
      y ~ 1 + a,
      data = data,
      family = family(),
      stanvars = family()$stanvars,
      backend = "cmdstanr",
      chains = num_parallel,
      cores = num_parallel,
      # seed = seed, #no seed for BRMS. BRMS RNG code is to be tested to, so would be not sensible
      silent = 2,
      refresh = 0,
      init = 0.1
    )
  }

  return(posterior_fit)
}

#' Check, that data of the posterior is close enough to the reference data.
#'
#' @param posterior_data Data fitted and drawn, a BRMS object, which is a list in R terms
#' @param name Name of the variable to check, as single string
#' @param reference Reference value to check against, single real scalar
#' @param thresh real scalar or 2-length vector of quantile bounds.
#' For scalar constructs bound as [tresh, 1-thresh]
#' thresh has to be inside the Unit-Interval.
#'
#' @return Single boolean succeess, fail or error
#' @export
#'
#' @examples ba_in <- 0.5
#' # with seed=1337 returns true in my case, but it may fail on other machines
#' # especially if the RNG generator does not work the same. In any case, you may check visually against the plot!
#' fit1 <- construct_brms(1000, ba = ba_in, 1.0, 2.0, exp, bayesim::betaprime, bayesim::rbetaprime, 2, seed=1337)
#' result <- test_brms_quantile(fit1, "b_a", ba_in, 0.025)
#' print(result)
#' plot(fit1)
test_brms_quantile <- function(posterior_data, name, reference, thresh, debug = FALSE) {
  if (!is.list(posterior_data)) {
    stop("The posterior_data frame has to be BRMS data, which itself is in R of type list")
  }
  if (!isSingleString(name)) {
    stop("The variable name argument has to be a single string")
  }
  if (!isNum_len(reference)) {
    stop("The reference data has to be a single real scalar")
  }
  if (isTRUE(any(thresh <= 0 | thresh >= 1))) {
    stop("The threshold has to be in the Unit-Interval")
  }
  if (!isLogic_len(debug)) {
    stop("The argument debug has to be a single boolean")
  }
  if (isNum_len(thresh, 1)) {
    if (thresh > 0.5) {
      stop("Any threshold scalar bigger 0.5 would create a bigger lower bound than upper bound!")
    }
    bounds <- c(thresh, 1 - thresh)
  } else if (isNum_len(thresh, 2)) {
    if(isFALSE(thresh[1] <= tresh[2])) {
      stop("If a 2 entry vector is used for the bounds, the first entry is the lower bound")
    }
    bounds <- thresh
  } else {
    stop("The quantile-thresholds can only be a scalar, or a 2 length vector")
  }

  calculated <- tryCatch({
     posterior::extract_variable_matrix(posterior_data, variable = name)
  }, error = function(e) {
    # if the extract fails, most probably cause is a wrong variable name string
    # instead of giving an unreadable error, throw clear warning and return false
    warning(paste0("In test_brms_quantile, the variable extraction of ", name, " failed.
                   Most probable cause, the variable-string was written wrong or did not exist somehow.
                   Return FALSE in this case."))
    warning(paste0("Original error was: ", e))
    return(FALSE)
  })

  quantiles <- unname(quantile(calculated, probs = bounds))
  if (debug) {
    print(paste("At: ", bounds, " the quantile is: ", quantiles))
    print(paste("The supplied reference value was: ", reference))
  }

  # if(quantiles[1] < reference && reference < quantiles[2])
  #   succeed()
  # else
  #   fail("The reference value was not within the quantiles of the given posterior data!")
  value <- quantiles[1] < reference && reference < quantiles[2]
  return(isTRUE(value))
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
#' @param save_folder The folder where to save the results. Default is set to
#' "tests/testthat/precalc_values/" for internal use.
#' @param x_int_prelink Interval of x-values, before the x_link function is used on them.
#' @param x_link Link function to be used on x_int_prelink. Default is identity.
#'
#' @details This function will generate (n_param^2)*n_x number of reference samples.
#' So be aware not to increase both too much. Otherwise you'll require alot of memory.
#' As an example, for tests in this repo using modest n_x=1000 and leaving n_param=10 as default,
#' this function already generates about 800KiB of files per PDF to be tested.
#'
#' @return Nothing, but saves input and reference files for later use!
#'
#' @examples eps <- 1e-6
#' bayesim:::density_lookup_generator(mu_int = c(eps, 1 - eps), shape_int = c(2, 10), x_int_prelink = c(eps, 1 - eps),
#'   density_fun = bayesim::dcauchitnormal, density_name = "cauchitnormal_demodata",
#'   save_folder = "")
#' # saves cauchit_normal_demodata_refdata and *_refpdf with density input and lookup-data respectively
#' readRDS("cauchitnormal_demodata_refdata")
density_lookup_generator <- function(n_param = 10, n_x = 100, eps = 10^-6, mu_int, shape_int, x_int_prelink,
                                     x_link = identity, density_fun, density_name, save_folder = "tests/testthat/precalc_values/") {

  # assert parameters (why all the assertions, if only I ever get to use it?)
  # dunno, just felt like it :P
  if (isFALSE(is.function(density_fun))) {
    stop("The given density_fun was no function. Expect density-function as input.")
  }
  if (!isSingleString(density_name)) {
    stop("The given density_name has to be a single string.")
  }
  if (!isSingleString(save_folder)) {
    stop("The given save_folder has to be a single string.")
  }
  if (!is.function(x_link)) {
    stop("The given x_link was no function. Expect link-function as input.")
  }
  if (!isNum_len(mu_int, 2)) {
    stop("The mu_int interval has to be a two long vector")
  }
  if (!isNum_len(shape_int, 2)) {
    stop("The shape_int interval has to be a two long vector")
  }
  if (!isNum_len(x_int_prelink, 2)) {
    stop("The x_int_prelink interval has to be a two long vector")
  }
  if (isFALSE(isInt_len(n_param) && n_param > 0)) {
    stop("The n_param has to be a single whole integer bigger 0")
  }
  if (isFALSE(isInt_len(n_x) && n_x > 0)) {
    stop("The n_x has to be a single whole integer bigger 0")
  }
  if (isFALSE(isNum_len(eps) && eps > 0)) {
    stop("The eps has to be a scalar bigger 0")
  }

  # generate input data
  mus <- seq(from = mu_int[1], to = mu_int[2], length.out = n_param)
  shapes <- seq(from = shape_int[1], to = shape_int[2], length.out = n_param)

  x_ref <- x_link(seq(from = x_int_prelink[1], to = x_int_prelink[2], length.out = n_x))
  y_ref <- matrix(data = list(), nrow = n_param, ncol = n_param)

  # generate reference data
  for (outer in 1:n_param) {
    for (inner in 1:n_param) {
      mu <- mus[outer]
      shape <- shapes[inner]
      y_ref[[outer, inner]] <- density_fun(x_ref, mu, shape)
    }
  }

  # pack data into frames
  inp_scalars <- data.frame(n = n_x, n_small = n_param, eps = eps)
  inp_vectors <- data.frame(mus = mus, shapes = shapes)
  x_data <- data.frame(x = x_ref)
  y_data <- data.frame(y = y_ref)
  return_data <- c(inp_scalars, inp_vectors, x_data)
  # then save them to RDS objects
  saveRDS(return_data, paste0(save_folder, density_name, "_refdata"))
  saveRDS(y_data, paste0(save_folder, density_name, "_refpdf"))
  # after reading the help, the issue with saveRDS (and save for that matter) is,
  # that it might be used incorrectly, which will not occur, if any bayesim dev implements them


  # If no errors occured, give feedback, of what density reference was saved where
  print(paste0(
    "Generated lookup data \"", density_name, "\" and saved to file in \"",
    "package-root/", save_folder, "\" folder for tests"
  ), quote = FALSE)
}
