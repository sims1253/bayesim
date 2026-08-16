# Model-bank default-prior guard -------------------------------------------------
#
# These tests run in the FAST tier (no brms/cmdstanr needed): the guard inside
# build_model_bank() fires BEFORE any brms call for the spec (template-data
# generation and grid parsing are brms-free; the first brms touch would be
# model_spec_hash() -> resolve_formula(), which sits AFTER the guard), so the
# error tests need no skip guard. BrmsFitter is a plain S7 class: constructing
# it does not require brms (verified empirically; it only touches brms inside
# fit_model()). The file-level brms skip below therefore comes AFTER this block
# so these tests always execute.

.describe_guard_generator <- function(data_spec, task_ctx) {
  n <- as.integer(data_spec$n %||% 10L)
  x <- stats::rnorm(n)
  y <- 2 * x + stats::rnorm(n)
  list(train = data.frame(y = y, x = x), response = "y")
}

.describe_guard_grid <- function() {
  grid <- data.frame(model = "gaussian", stringsAsFactors = FALSE)
  grid$formula <- list(y ~ x)
  grid$family <- list(NULL)
  grid
}

describe("model bank prior guard (build_model_bank)", {
  it("aborts with a bayesim_config_error when a spec has no explicit prior", {
    fitter <- BrmsFitter(precompile = TRUE) # allow_default_priors defaults FALSE

    expect_error(
      bayesim:::build_model_bank(
        fitter = fitter,
        fit_grid = .describe_guard_grid(),
        data_generator = .describe_guard_generator,
        data_spec_template = list(n = 10)
      ),
      class = "bayesim_config_error"
    )
  })

  it("names the offending spec (row index + deparsed formula) and the remedies", {
    fitter <- BrmsFitter(precompile = TRUE)

    err <- tryCatch(
      bayesim:::build_model_bank(
        fitter = fitter,
        fit_grid = .describe_guard_grid(),
        data_generator = .describe_guard_generator,
        data_spec_template = list(n = 10)
      ),
      error = function(e) e
    )
    expect_s3_class(err, "bayesim_config_error")

    msg <- conditionMessage(err)
    # Offending spec identified: fit_grid row index and deparsed formula.
    expect_match(msg, "fit_grid row 1", fixed = TRUE)
    expect_match(msg, "y ~ x", fixed = TRUE)
    # All three remedies are listed.
    expect_match(msg, "explicit", fixed = TRUE)
    expect_match(msg, "prior", fixed = TRUE)
    expect_match(msg, "precompile = FALSE", fixed = TRUE)
    expect_match(msg, "allow_default_priors = TRUE", fixed = TRUE)
    # The hazard itself is explained.
    expect_match(msg, "template", ignore.case = TRUE)
  })

  it("passes the prior guard when allow_default_priors = TRUE", {
    bayesim:::.reset_warn_once()
    fitter <- BrmsFitter(
      precompile = TRUE,
      allow_default_priors = TRUE
    )

    # The opt-in notice must fire exactly once as a warning...
    noticed <- FALSE
    outcome <- tryCatch(
      withCallingHandlers(
        bayesim:::build_model_bank(
          fitter = fitter,
          fit_grid = .describe_guard_grid(),
          data_generator = .describe_guard_generator,
          data_spec_template = list(n = 10)
        ),
        warning = function(w) {
          if (
            grepl("default priors", conditionMessage(w), ignore.case = TRUE)
          ) {
            noticed <<- TRUE
          }
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )

    # ...and the guard must NOT have fired. Beyond the guard, the outcome is
    # environment-dependent: with brms + a Stan toolchain the bank builds; on a
    # machine without brms, model_spec_hash() -> resolve_formula() fails with a
    # plain "no package called 'brms'" error; with brms but no toolchain the
    # compile failure surfaces as a bayesim_internal_error. All of those prove
    # the guard passed; only the guard's own bayesim_config_error would not.
    expect_true(noticed)
    if (inherits(outcome, "error")) {
      expect_false(
        grepl("no explicit prior", conditionMessage(outcome), fixed = TRUE),
        info = paste0(
          "prior guard fired despite allow_default_priors = TRUE: ",
          conditionMessage(outcome)
        )
      )
      expect_false(
        inherits(outcome, "bayesim_config_error"),
        info = paste0(
          "unexpected config error past the prior guard: ",
          conditionMessage(outcome)
        )
      )
    } else {
      expect_type(outcome, "list")
    }
  })
})

# D3: brms_model() / model_grid() ergonomics.
skip_if_not(requireNamespace("brms", quietly = TRUE), "brms not available")

describe("brms_model", {
  it("builds a one-row spec from a formula and family", {
    spec <- brms_model(y ~ x, family = gaussian(), prior = NULL)
    expect_type(spec, "list")
    expect_equal(
      names(spec),
      c("formula", "family", "prior", "stanvars")
    )
    expect_equal(spec$formula, y ~ x)
    expect_false(is.null(spec$family))
  })

  it("rejects a non-formula", {
    expect_error(brms_model("not a formula"), class = "bayesim_config_error")
  })
})

describe("model_grid", {
  it("assembles a tibble of named specs", {
    grid <- model_grid(
      gaussian = brms_model(y ~ x, gaussian()),
      student = brms_model(y ~ x, brms::brmsfamily("student"))
    )
    expect_s3_class(grid, "data.frame")
    expect_equal(nrow(grid), 2L)
    expect_equal(grid$model, c("gaussian", "student"))
    expect_true(all(
      c("formula", "family", "prior", "stanvars") %in% names(grid)
    ))
    # formula list-column carries the specs.
    expect_equal(grid$formula[[1]], y ~ x)
    # family list-column is populated.
    expect_false(is.null(grid$family[[2]]))
  })

  it("errors without names", {
    expect_error(
      model_grid(brms_model(y ~ x, gaussian())),
      class = "bayesim_config_error"
    )
  })

  it("errors on a non-brms_model spec", {
    expect_error(
      model_grid(bad = list(formula = y ~ x)),
      class = "bayesim_config_error"
    )
  })
})
