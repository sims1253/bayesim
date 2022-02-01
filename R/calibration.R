#' Title
#'
#' @param model_name
#' @param results
#' @param samples
#'
#' @return
#' @export
#'
#' @examples
calibration_result_text <- function(model_name, results, samples) {
  print(model_name)
  greater_lower_ci <- unlist(
    lapply(
      results$greater,
      function(x) {
        x$hypothesis$CI.Lower
      }
    )
  )
  print(
    paste(
      length(greater_lower_ci[greater_lower_ci < 0]),
      "/",
      samples,
      " of 'x>0' lower CI bounds are < 0",
      sep = ""
    )
  )

  equal_lower_ci <- unlist(
    lapply(
      results$equal,
      function(x) {
        x$hypothesis$CI.Lower
      }
    )
  )
  print(
    paste(
      length(equal_lower_ci[equal_lower_ci < 0]),
      "/",
      samples,
      " of 'x=0' lower CI bounds are < 0",
      sep = ""
    )
  )
  print("=================================================")
}

#' Title
#'
#' @param dataset
#' @param prefit
#'
#' @return
#'
#' @examples
parallel_run <- function(dataset, prefit) {
  fit <- stats::update(prefit,
    newdata = dataset,
    refresh = 0,
    silent = 2,
    chains = 4,
    control = list(adapt_delta = 0.9),
    prior = c(
      brms::prior("", class = "Intercept")
    )
  )
  return(
    list(
      greater = brms::hypothesis(fit, hypothesis = "x>0"),
      equal = brms::hypothesis(fit, hypothesis = "x=0")
    )
  )
}

#' Title
#'
#' @param data_list
#' @param brms_formular
#' @param ncores
#' @param brms_family
#' @param stanvars
#' @param exports
#'
#' @return
#' @export
#'
#' @examples
calibrate_formula <- function(data_list,
                              brms_formular,
                              ncores,
                              brms_family,
                              stanvars,
                              exports) {
  prefit <- brms::brm(
    brms_formular,
    data = data_list[[1]],
    family = brms_family,
    stanvars = stanvars,
    chains = 0,
    refresh = 0,
    silent = 2,
    control = list(adapt_delta = 0.9),
    prior = c(
      brms::prior("", class = "Intercept")
    )
  )

  cluster <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cluster)
  `%dopar%` <- foreach::`%dopar%`

  results <- foreach::foreach(
    dataset = data_list,
    .packages = c("brms"),
    .export = exports
  ) %dopar% {
    parallel_run(dataset, prefit)
  }
  parallel::stopCluster(cluster)

  return(
    list(
      greater = lapply(results, function(x) {
        x$greater
      }),
      equal = lapply(results, function(x) {
        x$equal
      })
    )
  )
}

#' Title
#'
#' @param dataset_N
#' @param ncores
#' @param brms_family
#' @param stanvars
#' @param exports
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
calibrate_family <- function(dataset_N,
                             ncores,
                             brms_family,
                             stanvars,
                             exports,
                             ...) {
  data_list <- vector("list", length = dataset_N)
  for (i in seq_len(dataset_N)) {
    data_list[i] <- list(basedag_data(...))
  }
  hist(data_list[[1]]$y, main = "data", xlab = "data")
  true_model_logit_sim <- calibrate_formula(data_list,
    brms_formular = brms::brmsformula(y ~ x + z1 + z2),
    ncores = ncores,
    brms_family = brms_family,
    stanvars = stanvars,
    exports = exports
  )
  calibration_result_text(model_name = "True Model", results = true_model_logit_sim, samples = dataset_N)

  z1_model_logit_sim <- calibrate_formula(data_list,
    brms_formular = brms::brmsformula(y ~ x + z2),
    ncores = ncores,
    brms_family = brms_family,
    stanvars = stanvars,
    exports = exports
  )
  calibration_result_text(model_name = "z1 Model", results = z1_model_logit_sim, samples = dataset_N)

  z2_model_logit_sim <- calibrate_formula(data_list,
    brms_formular = brms::brmsformula(y ~ x + z1),
    ncores = ncores,
    brms_family = brms_family,
    stanvars = stanvars,
    exports = exports
  )
  calibration_result_text(model_name = "z2 Model", results = z2_model_logit_sim, samples = dataset_N)

  z3_model_logit_sim <- calibrate_formula(data_list,
    brms_formular = brms::brmsformula(y ~ x + z1 + z2 + z3),
    ncores = ncores,
    brms_family = brms_family,
    stanvars = stanvars,
    exports = exports
  )
  calibration_result_text(model_name = "z3 Model", results = z3_model_logit_sim, samples = dataset_N)

  z4_model_logit_sim <- calibrate_formula(data_list,
    brms_formular = brms::brmsformula(y ~ x + z1 + z2 + z4),
    ncores = ncores,
    brms_family = brms_family,
    stanvars = stanvars,
    exports = exports
  )
  calibration_result_text(model_name = "z4 Model", results = z4_model_logit_sim, samples = dataset_N)
}
