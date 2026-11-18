# Session model cache for CmdStanFitter. The compiler is mocked, so these
# tests need cmdstanr installed but never compile a Stan program.
skip_if_not_installed("cmdstanr")

test_that("editing a Stan file invalidates the session model cache", {
  withr::local_options(list(bayesim.cmdstan_models = NULL))
  path <- withr::local_tempfile(fileext = ".stan")
  writeLines("parameters { real theta; } model { theta ~ normal(0, 1); }", path)
  compiled <- character()
  local_mocked_bindings(
    cmdstan_model = function(stan_file = NULL, ...) {
      compiled <<- c(compiled, paste(readLines(stan_file), collapse = "\n"))
      list(code = compiled[[length(compiled)]])
    },
    .package = "cmdstanr"
  )
  fitter <- CmdStanFitter(
    stan_file = path,
    stan_data = function(data_bundle, fit_spec) list()
  )

  first <- bayesim:::.cmdstan_model(fitter)
  # Identical contents are a cache hit: no second compile.
  again <- bayesim:::.cmdstan_model(fitter)
  expect_identical(again, first)
  expect_length(compiled, 1L)

  # Changed contents at the same path must compile the new program.
  writeLines("parameters { real theta; } model { theta ~ normal(10, 1); }", path)
  second <- bayesim:::.cmdstan_model(fitter)
  expect_length(compiled, 2L)
  expect_false(identical(first$code, second$code))
  expect_match(second$code, "normal(10, 1)", fixed = TRUE)
})
