library(bayesim)

# Package-external metric extension interface: direct subclass construction
# (Metric is abstract, so subclass constructors honor declared defaults) and
# conformance validation (validate_metric() representative execution path).

describe("direct subclass construction", {
  it("rejects direct instantiation of the abstract Metric base", {
    expect_error(Metric(), "abstract")
  })

  it("honors subclass-declared defaults for inherited Metric properties", {
    # Metric is abstract, so S7's generated subclass constructor uses the
    # subclass's own redeclarations of name/needs/required/schema instead of
    # silently delegating to the parent's defaults. External metrics can be
    # constructed directly: MyMetric() retains what MyMetric declared.
    ExtDefaultsMetric <- S7::new_class(
      "ExternalDefaultsMetric",
      parent = Metric,
      properties = list(
        name = S7::new_property(
          S7::class_character,
          default = "ext_defaults"
        ),
        needs = S7::new_property(
          S7::class_character,
          default = c("predictions", "log_lik")
        ),
        required = S7::new_property(S7::class_logical, default = TRUE),
        schema = S7::new_property(
          S7::class_list,
          default = list(value = list(role = "estimate"))
        )
      )
    )

    metric <- ExtDefaultsMetric()

    expect_true(S7::S7_inherits(metric, Metric))
    # S7 qualifies class names created inside a namespace (e.g.
    # "bayesim::ExternalDefaultsMetric" under test_local), so assert the
    # concrete subclass via the class object rather than its name string.
    expect_true(S7::S7_inherits(metric, ExtDefaultsMetric))
    expect_equal(metric@name, "ext_defaults")
    expect_equal(metric@needs, c("predictions", "log_lik"))
    expect_true(metric@required)
    expect_equal(metric@summary_type, "mean")
    expect_equal(metric@schema, list(value = list(role = "estimate")))
  })

  it("lets explicit arguments override subclass defaults", {
    ExtDefaultsMetric <- S7::new_class(
      "ExternalOverrideMetric",
      parent = Metric,
      properties = list(
        name = S7::new_property(
          S7::class_character,
          default = "ext_override"
        ),
        needs = S7::new_property(
          S7::class_character,
          default = "predictions"
        )
      )
    )

    metric <- ExtDefaultsMetric(name = "custom_name")

    expect_equal(metric@name, "custom_name")
    expect_equal(metric@needs, "predictions")
  })

  it("falls back to Metric defaults when nothing is declared", {
    PlainMetric <- S7::new_class("ExternalPlainMetric", parent = Metric)

    metric <- PlainMetric(name = "plain")

    expect_equal(metric@name, "plain")
    expect_equal(metric@needs, character())
    expect_false(metric@required)
    expect_equal(metric@summary_type, "mean")
    expect_equal(metric@schema, list())
  })

  it("forwards subclass-specific properties", {
    ThresholdMetric <- S7::new_class(
      "ExternalThresholdMetric",
      parent = Metric,
      properties = list(
        threshold = S7::new_property(S7::class_double, default = 0.5)
      )
    )

    expect_equal(ThresholdMetric(name = "thr")@threshold, 0.5)
    expect_equal(
      ThresholdMetric(name = "thr", threshold = 0.9)@threshold,
      0.9
    )
  })

  it("validates declarations via validate_metric()", {
    PlainMetric <- S7::new_class("ExternalUnnamedMetric", parent = Metric)

    # No name supplied and no subclass-declared default: the name stays
    # empty, which validate_metric() rejects.
    expect_error(validate_metric(PlainMetric()), "name")
    expect_error(
      validate_metric(PlainMetric(name = "has__separator")),
      "__"
    )
  })
})

describe("validate_metric() representative execution path", {
  RepresentativeMetric <- S7::new_class(
    "ExternalRepresentativeMetric",
    parent = Metric,
    properties = list(
      name = S7::new_property(
        S7::class_character,
        default = "ext_representative"
      ),
      needs = S7::new_property(
        S7::class_character,
        default = "predictions"
      )
    )
  )
  S7::method(compute_metric, RepresentativeMetric) <- function(
    metric,
    fit_result,
    data_bundle,
    context,
    task_ctx
  ) {
    list(
      value = mean(context$predictions$predicted_mean),
      n_obs = length(context$predictions$predicted_mean)
    )
  }

  representative_fixtures <- function() {
    n_test <- 4L
    samples <- matrix(
      stats::rnorm(8L * n_test),
      nrow = 8L,
      ncol = n_test
    )
    list(
      fit_result = new_fit_result(
        success = TRUE,
        fit = list(),
        draws = cbind(intercept = rep(1, 8L), slope = rep(2, 8L)),
        diagnostics = list(),
        timing = list(total = 0)
      ),
      data_bundle = list(
        train = data.frame(
          outcome = seq(1, 2, length.out = 8L),
          feature = seq(-1, 1, length.out = 8L)
        ),
        test = data.frame(
          outcome = seq(1, 2, length.out = n_test),
          feature = seq(-0.5, 0.5, length.out = n_test)
        ),
        response = "outcome",
        true_params = c(intercept = 1, slope = 2, sigma = 1),
        vars_of_interest = c("intercept", "slope", "sigma"),
        meta = list()
      ),
      context = list(
        predictions = list(
          predicted_mean = rep(1.5, n_test),
          predicted_samples = samples,
          predicted_sd = apply(samples, 2, stats::sd)
        )
      ),
      task_ctx = list(
        task_id = "validate_metric_test",
        data_idx = 1L,
        fit_idx = 1L,
        rep_idx = 1L,
        seed = 42L
      )
    )
  }

  it("executes compute_metric with caller-supplied representative values", {
    fx <- representative_fixtures()
    metric <- RepresentativeMetric()

    expect_silent(validate_metric(
      metric,
      fit_result = fx$fit_result,
      data_bundle = fx$data_bundle,
      context = fx$context,
      task_ctx = fx$task_ctx
    ))
    expect_equal(
      validate_metric(
        metric,
        fit_result = fx$fit_result,
        data_bundle = fx$data_bundle,
        context = fx$context,
        task_ctx = fx$task_ctx
      )@name,
      "ext_representative"
    )
  })

  it("supplies defaults for context and task_ctx", {
    fx <- representative_fixtures()

    # context and task_ctx are optional; the metric degrades to NA without
    # crashing the validator's plumbing itself.
    DegradingMetric <- S7::new_class(
      "ExternalDegradingMetric",
      parent = Metric,
      properties = list(
        name = S7::new_property(
          S7::class_character,
          default = "ext_degrading"
        )
      )
    )
    S7::method(compute_metric, DegradingMetric) <- function(
      metric,
      fit_result,
      data_bundle,
      context,
      task_ctx
    ) {
      list(value = if (is.null(context$predictions)) NA_real_ else 1)
    }

    expect_silent(validate_metric(
      DegradingMetric(),
      fit_result = fx$fit_result,
      data_bundle = fx$data_bundle
    ))
  })

  it("rejects malformed compute_metric output", {
    BadOutputMetric <- S7::new_class(
      "ExternalBadOutputMetric",
      parent = Metric,
      properties = list(
        name = S7::new_property(
          S7::class_character,
          default = "ext_bad_output"
        )
      )
    )
    S7::method(compute_metric, BadOutputMetric) <- function(
      metric,
      fit_result,
      data_bundle,
      context,
      task_ctx
    ) {
      list(nested = list(a = 1))
    }

    fx <- representative_fixtures()
    expect_error(
      validate_metric(
        BadOutputMetric(),
        fit_result = fx$fit_result,
        data_bundle = fx$data_bundle
      ),
      class = "bayesim_validation_error"
    )
  })

  it("flags compute_metric failures as validation errors", {
    FailingMetric <- S7::new_class(
      "ExternalFailingMetric",
      parent = Metric,
      properties = list(
        name = S7::new_property(
          S7::class_character,
          default = "ext_failing"
        )
      )
    )
    S7::method(compute_metric, FailingMetric) <- function(
      metric,
      fit_result,
      data_bundle,
      context,
      task_ctx
    ) {
      stop("deliberate representative failure")
    }

    fx <- representative_fixtures()
    expect_error(
      validate_metric(
        FailingMetric(),
        fit_result = fx$fit_result,
        data_bundle = fx$data_bundle
      ),
      "deliberate representative failure",
      class = "bayesim_validation_error"
    )
  })

  it("checks declared schema fields against produced output", {
    SchemaMetric <- S7::new_class(
      "ExternalSchemaMetric",
      parent = Metric,
      properties = list(
        name = S7::new_property(
          S7::class_character,
          default = "ext_schema"
        ),
        schema = S7::new_property(
          S7::class_list,
          default = list(
            value = list(role = "estimate"),
            ghost_field = list(role = "estimate")
          )
        )
      )
    )
    S7::method(compute_metric, SchemaMetric) <- function(
      metric,
      fit_result,
      data_bundle,
      context,
      task_ctx
    ) {
      list(value = 1)
    }

    fx <- representative_fixtures()
    expect_error(
      validate_metric(
        SchemaMetric(),
        fit_result = fx$fit_result,
        data_bundle = fx$data_bundle
      ),
      "ghost_field",
      class = "bayesim_validation_error"
    )
  })

  it("stays declaration-only when no representative values are supplied", {
    NoComputeMetric <- S7::new_class(
      "ExternalNoComputeMetric",
      parent = Metric,
      properties = list(
        name = S7::new_property(
          S7::class_character,
          default = "ext_no_compute"
        )
      )
    )
    # No compute_metric method is registered; declaration-only validation
    # must keep passing (compute is only exercised when representative values
    # are supplied).
    expect_silent(validate_metric(NoComputeMetric()))
  })

  it("requires fit_result and data_bundle together", {
    fx <- representative_fixtures()
    metric <- RepresentativeMetric()

    expect_error(
      validate_metric(metric, fit_result = fx$fit_result),
      "data_bundle"
    )
    expect_error(
      validate_metric(metric, data_bundle = fx$data_bundle),
      "fit_result"
    )
  })

  it("requires a bayesim_fit_result representative", {
    fx <- representative_fixtures()
    metric <- RepresentativeMetric()

    expect_error(
      validate_metric(
        metric,
        fit_result = list(not_a = "fit_result"),
        data_bundle = fx$data_bundle,
        context = fx$context
      ),
      "bayesim_fit_result"
    )
  })
})
