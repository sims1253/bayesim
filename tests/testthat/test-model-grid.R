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
