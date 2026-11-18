library(bayesim)

test_that("formula list-columns participate in the public config fingerprint", {
  generator <- function(data_spec, task_ctx) {
    list(
      train = data.frame(y = 1, x = 1),
      test = NULL,
      response = "y",
      true_params = c(x = 1),
      vars_of_interest = "x"
    )
  }
  make_config <- function(formula) {
    fit_grid <- data.frame(model = "gaussian")
    fit_grid$formula <- list(formula)
    simulation_config(
      data_grid = data.frame(n = 1L),
      fit_grid = fit_grid,
      data_generator = generator,
      fitter = LinearRegressionFitter(n_draws = 10L),
      metrics = list(),
      n_replicates = 1L,
      seed = 101L
    )
  }

  expect_false(identical(
    config_fingerprint(make_config(y ~ x)),
    config_fingerprint(make_config(y ~ 1))
  ))
})

test_that("generator default arguments participate in the config fingerprint", {
  gen1 <- function(data_spec, task_ctx, sigma = 1) {
    x <- stats::rnorm(data_spec$n)
    y <- x + stats::rnorm(data_spec$n, sd = sigma)
    list(train = data.frame(y = y, x = x), response = "y")
  }
  gen2 <- gen1
  formals(gen2)$sigma <- 10 # same body, same arg names, different default
  gen3 <- function(data_spec, task_ctx, sigma = 1) {
    x <- stats::rnorm(data_spec$n)
    y <- x + stats::rnorm(data_spec$n, sd = sigma)
    list(train = data.frame(y = y, x = x), response = "y")
  }
  gen4 <- local({
    mu <- 1
    function(data_spec, task_ctx, center = mu) stats::rnorm(data_spec$n, center)
  })
  gen5 <- local({
    mu <- 5
    function(data_spec, task_ctx, center = mu) stats::rnorm(data_spec$n, center)
  })
  make_config <- function(generator) {
    simulation_config(
      data_grid = data.frame(n = 20L),
      fit_grid = data.frame(model = "lm"),
      data_generator = generator,
      fitter = LinearRegressionFitter(n_draws = 10L),
      metrics = list(),
      n_replicates = 1L,
      seed = 42L
    )
  }

  gen_body <- function(data_spec, task_ctx, sigma = 1) {
    x <- stats::rnorm(data_spec$n)
    y <- 2 * x + stats::rnorm(data_spec$n, sd = sigma)
    list(train = data.frame(y = y, x = x), response = "y")
  }

  # A different body must not share a fingerprint...
  expect_false(identical(
    config_fingerprint(make_config(gen1)),
    config_fingerprint(make_config(gen_body))
  ))
  # ...nor may differing default arguments...
  expect_false(identical(
    config_fingerprint(make_config(gen1)),
    config_fingerprint(make_config(gen2))
  ))
  # ...and a closure binding used only as a default value must participate
  # like a body free variable.
  expect_false(identical(
    config_fingerprint(make_config(gen4)),
    config_fingerprint(make_config(gen5))
  ))
  # The same definition yields the same fingerprint across calls.
  expect_identical(
    config_fingerprint(make_config(gen1)),
    config_fingerprint(make_config(gen3))
  )
})

test_that("checkpoint recovery falls back after checksum corruption", {
  result_path <- file.path(withr::local_tempdir(), "checksum-fallback")
  fingerprint <- "focused-corruption-test"
  init_checkpoint_dir(result_path, config_fingerprint = fingerprint)

  first <- new_task_result(
    task_id = "d001_f001_r00001",
    status = "success",
    metrics = list(value = 1)
  )
  first_grid <- data.frame(
    task_id = first$task_id,
    status = first$status,
    stop_reason = NA_character_
  )
  write_checkpoint(
    result_path,
    task_grid = first_grid,
    task_results = list(first),
    config_fingerprint = fingerprint
  )

  second <- new_task_result(
    task_id = "d001_f001_r00002",
    status = "success",
    metrics = list(value = 2)
  )
  second_grid <- rbind(
    first_grid,
    data.frame(
      task_id = second$task_id,
      status = second$status,
      stop_reason = NA_character_
    )
  )
  write_checkpoint(
    result_path,
    task_grid = second_grid,
    task_results = list(second),
    config_fingerprint = fingerprint
  )

  latest_meta <- file.path(
    result_path,
    "checkpoints",
    "cp_000002",
    "meta.json"
  )
  expect_true(file.remove(latest_meta))

  recovered <- NULL
  expect_warning(
    recovered <- get_latest_valid_checkpoint(
      result_path,
      config_fingerprint = fingerprint
    ),
    "invalid checksums"
  )
  expect_equal(recovered$checkpoint_id, 1L)
  expect_equal(recovered$results_df$value, 1)
})
