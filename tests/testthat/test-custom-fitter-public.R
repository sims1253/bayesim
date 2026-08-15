library(bayesim)

describe("package-external extension workflow", {
  it("composes a generator, fitter, and metric through run_simulation", {
    # Every bayesim symbol used by this extension is part of the public API.
    extension_api <- c(
      "Fitter",
      "Metric",
      "fit_model",
      "extract_draws",
      "predict_fit",
      "log_lik_matrix",
      "fit_diagnostics",
      "compute_metric",
      "new_fit_result",
      "validate_fitter",
      "simulation_config",
      "run_simulation"
    )
    expect_setequal(
      intersect(extension_api, getNamespaceExports("bayesim")),
      extension_api
    )

    ExternalGaussianFitter <- S7::new_class(
      "ExternalGaussianFitter",
      parent = Fitter,
      properties = list(
        name = S7::new_property(
          S7::class_character,
          default = "external_gaussian"
        ),
        supports_predictions = S7::new_property(
          S7::class_logical,
          default = TRUE
        ),
        supports_log_lik = S7::new_property(
          S7::class_logical,
          default = TRUE
        )
      )
    )

    S7::method(fit_model, ExternalGaussianFitter) <- function(
      fitter,
      data_bundle,
      fit_spec,
      seed,
      task_ctx
    ) {
      fit <- stats::lm(outcome ~ feature, data = data_bundle$train)
      center <- stats::coef(fit)
      offsets <- seq(-0.1, 0.1, length.out = 11L)
      draws <- cbind(
        intercept = center[[1L]] + offsets,
        slope = center[[2L]] + offsets,
        sigma = rep(1, length(offsets))
      )
      new_fit_result(
        fit = list(model = fit, data_bundle = data_bundle),
        draws = draws,
        diagnostics = list(converged = TRUE),
        timing = list(total = 0),
        data_bundle = data_bundle
      )
    }

    S7::method(extract_draws, ExternalGaussianFitter) <- function(
      fitter,
      fit_result,
      variables = NULL
    ) {
      draws <- fit_result$draws
      if (!is.null(variables)) {
        draws <- draws[, intersect(variables, colnames(draws)), drop = FALSE]
      }
      draws
    }

    S7::method(predict_fit, ExternalGaussianFitter) <- function(
      fitter,
      fit_result,
      newdata = NULL,
      seed = NULL
    ) {
      data <- if (is.null(newdata)) {
        fit_result$data_bundle$train
      } else {
        newdata
      }
      draws <- fit_result$draws
      samples <- outer(draws[, "intercept"], rep(1, nrow(data))) +
        outer(draws[, "slope"], data$feature)
      list(
        predicted_mean = colMeans(samples),
        predicted_samples = samples,
        predicted_sd = apply(samples, 2, stats::sd)
      )
    }

    S7::method(log_lik_matrix, ExternalGaussianFitter) <- function(
      fitter,
      fit_result,
      newdata = NULL
    ) {
      data <- if (is.null(newdata)) {
        fit_result$data_bundle$train
      } else {
        newdata
      }
      predictions <- predict_fit(fitter, fit_result, newdata = data)
      residuals <- sweep(
        predictions$predicted_samples,
        2L,
        data$outcome,
        "-"
      )
      -0.5 * log(2 * pi) - 0.5 * residuals^2
    }

    ExternalContractMetric <- S7::new_class(
      "ExternalContractMetric",
      parent = Metric,
      properties = list(
        name = S7::new_property(
          S7::class_character,
          default = "external_contract"
        ),
        needs = S7::new_property(
          S7::class_character,
          default = c("predictions", "log_lik")
        ),
        required = S7::new_property(S7::class_logical, default = TRUE)
      )
    )

    S7::method(compute_metric, ExternalContractMetric) <- function(
      metric,
      fit_result,
      data_bundle,
      context,
      task_ctx
    ) {
      list(
        prediction_mean = mean(context$predictions$predicted_mean),
        mean_log_lik = mean(context$log_lik),
        prediction_n = ncol(context$predictions$predicted_samples),
        log_lik_n = ncol(context$log_lik)
      )
    }

    external_generator <- function(data_spec, task_ctx) {
      train_x <- seq(-1, 1, length.out = data_spec$n_train)
      test_x <- seq(-0.75, 0.75, length.out = data_spec$n_test)
      list(
        train = data.frame(
          outcome = 1 + 2 * train_x,
          feature = train_x
        ),
        test = data.frame(
          outcome = 1 + 2 * test_x,
          feature = test_x
        ),
        response = "outcome",
        true_params = c(intercept = 1, slope = 2, sigma = 1),
        vars_of_interest = c("intercept", "slope", "sigma"),
        meta = list(generator = "external")
      )
    }

    fitter <- ExternalGaussianFitter()
    representative_bundle <- external_generator(
      list(n_train = 12L, n_test = 5L),
      list(task_id = "conformance")
    )
    expect_silent(validate_fitter(
      fitter,
      smoke_test = TRUE,
      data_bundle = representative_bundle,
      fit_spec = list(model = "gaussian")
    ))

    config <- simulation_config(
      data_grid = data.frame(n_train = 12L, n_test = 5L),
      fit_grid = data.frame(model = "gaussian"),
      data_generator = external_generator,
      fitter = fitter,
      metrics = list(ExternalContractMetric(
        name = "external_contract",
        needs = c("predictions", "log_lik"),
        required = TRUE
      )),
      n_replicates = 2L,
      seed = 42L
    )
    result <- run_simulation(
      config,
      resume = "never",
      progress = FALSE,
      verbose = FALSE
    )

    expect_true(all(result$summary$status == "success"))
    expect_equal(
      result$summary$external_contract__prediction_n,
      rep(5, 2L)
    )
    expect_equal(
      result$summary$external_contract__log_lik_n,
      rep(5, 2L)
    )
    expect_true(all(is.finite(
      result$summary$external_contract__prediction_mean
    )))
    expect_true(all(is.finite(
      result$summary$external_contract__mean_log_lik
    )))
  })

  it("defaults optional fitter capabilities to unsupported", {
    CoreOnlyFitter <- S7::new_class(
      "ExternalCoreOnlyFitter",
      parent = Fitter,
      properties = list(
        name = S7::new_property(S7::class_character, default = "core_only")
      )
    )
    expect_false(CoreOnlyFitter()@supports_predictions)
    expect_false(CoreOnlyFitter()@supports_log_lik)
    expect_false(CoreOnlyFitter()@supports_loo)
    expect_equal(fit_diagnostics(CoreOnlyFitter(), new_fit_result()), list())
  })

  it("keeps external fitter method dispatch alive with workers = 1", {
    # Dogfood regression: an external S7 fitter that passes
    # validate_fitter(smoke_test = TRUE) used to lose fit_model() dispatch
    # when run_simulation(workers = 1) shipped tasks to a single isolated
    # mirai daemon — S7 method tables are registered per process, so the
    # daemon could not see the method defined in the controller session.
    # workers = 1 must be genuinely sequential; only workers >= 2 sets up
    # daemons (parallel semantics preserved).
    ExternalSequentialFitter <- S7::new_class(
      "ExternalSequentialFitter",
      parent = Fitter,
      properties = list(
        name = S7::new_property(
          S7::class_character,
          default = "external_sequential"
        )
      )
    )

    S7::method(fit_model, ExternalSequentialFitter) <- function(
      fitter,
      data_bundle,
      fit_spec,
      seed,
      task_ctx
    ) {
      draws <- cbind(
        intercept = rep(1, 4L),
        slope = rep(2, 4L),
        sigma = rep(1, 4L)
      )
      new_fit_result(
        fit = list(source = "external"),
        draws = draws,
        diagnostics = list(converged = TRUE),
        timing = list(total = 0),
        data_bundle = data_bundle
      )
    }

    S7::method(extract_draws, ExternalSequentialFitter) <- function(
      fitter,
      fit_result,
      variables = NULL
    ) {
      fit_result$draws
    }

    sequential_generator <- function(data_spec, task_ctx) {
      train_x <- seq(-1, 1, length.out = data_spec$n_train)
      test_x <- seq(-0.5, 0.5, length.out = data_spec$n_test)
      list(
        train = data.frame(
          outcome = 1 + 2 * train_x,
          feature = train_x
        ),
        test = data.frame(
          outcome = 1 + 2 * test_x,
          feature = test_x
        ),
        response = "outcome",
        true_params = c(intercept = 1, slope = 2, sigma = 1),
        vars_of_interest = c("intercept", "slope", "sigma"),
        meta = list(generator = "external")
      )
    }

    fitter <- ExternalSequentialFitter()
    representative_bundle <- sequential_generator(
      list(n_train = 8L, n_test = 4L),
      list(task_id = "conformance")
    )
    expect_silent(validate_fitter(
      fitter,
      smoke_test = TRUE,
      data_bundle = representative_bundle,
      fit_spec = list(model = "gaussian")
    ))

    config <- simulation_config(
      data_grid = data.frame(n_train = 8L, n_test = 4L),
      fit_grid = data.frame(model = "gaussian"),
      data_generator = sequential_generator,
      fitter = fitter,
      metrics = list(),
      n_replicates = 1L,
      seed = 42L
    )

    result <- run_simulation(
      config,
      resume = "never",
      progress = FALSE,
      verbose = FALSE,
      workers = 1
    )

    expect_false(mirai::daemons_set())
    expect_true(all(result$summary$status == "success"))
    expect_equal(nrow(result$errors), 0L)
  })
})
