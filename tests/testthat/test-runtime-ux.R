# F1-F4: runtime UX helpers.
.gen <- function(data_spec, task_ctx) {
  n <- data_spec$n
  x <- stats::rnorm(n)
  y <- x + stats::rnorm(n)
  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(x = 1),
    vars_of_interest = "x",
    meta = list()
  )
}

describe("F1 preflight", {
  it("reports task count, grid shape, and unmet needs", {
    config <- simulation_config(
      data_grid = data.frame(n = c(20, 40)),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 50L),
      metrics = list(pred_rmse_metric()),
      n_replicates = 2L,
      seed = 1L
    )
    info <- preflight(config, condensed = FALSE)
    expect_equal(info$n_tasks, 4L)
    expect_equal(info$data_grid, 2L)
    expect_equal(info$fit_grid, 1L)
    expect_equal(info$n_replicates, 2L)
  })

  it("condensed mode runs without error", {
    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 50L),
      metrics = list(),
      n_replicates = 1L,
      seed = 1L
    )
    expect_invisible(preflight(config, condensed = TRUE))
  })

  it("counts distinct brms model specs rather than fit-grid rows", {
    skip_if_not(requireNamespace("brms", quietly = TRUE))
    grid <- data.frame(model = c("a", "a-copy"))
    grid$formula <- list(y ~ x, y ~ x)
    grid$family <- list(stats::gaussian(), stats::gaussian())
    config <- simulation_config(
      data_grid = data.frame(n = 20L),
      fit_grid = grid,
      data_generator = .gen,
      fitter = BrmsFitter(),
      metrics = list(),
      n_replicates = 1L,
      seed = 1L
    )

    info <- suppressMessages(preflight(config, condensed = TRUE))

    expect_equal(info$n_compile, 1L)
  })
})

describe("F2 failed_tasks", {
  it("returns an empty tibble when all tasks succeed", {
    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 50L),
      metrics = list(),
      n_replicates = 2L,
      seed = 1L
    )
    res <- run_simulation(config, resume = "never", progress = FALSE)
    ft <- failed_tasks(res)
    expect_s3_class(ft, "data.frame")
    expect_equal(nrow(ft), 0L)
  })
})

describe("F3 as_tibble", {
  it("returns the summary tibble", {
    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 50L),
      metrics = list(posterior_summary_metric()),
      n_replicates = 2L,
      seed = 1L
    )
    res <- run_simulation(config, resume = "never", progress = FALSE)
    tb <- tibble::as_tibble(res)
    expect_s3_class(tb, "tbl_df")
    expect_equal(nrow(tb), 2L)
  })
})

describe("F4 seed error message", {
  it("errors with a helpful message when seed is missing", {
    expect_error(
      simulation_config(
        data_grid = data.frame(n = 20),
        fit_grid = data.frame(model = "lm"),
        data_generator = .gen,
        fitter = LinearRegressionFitter(n_draws = 50L)
      ),
      "RNG stream"
    )
  })
})

describe("F5 run reporting", {
  it("separates verbosity from progress and prints resumable paths", {
    path <- withr::local_tempdir()
    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 20L),
      metrics = list(),
      n_replicates = 1L,
      seed = 11L,
      result_path = file.path(path, "run")
    )
    result <- expect_silent(
      run_simulation(
        config,
        resume = "never",
        progress = FALSE,
        verbose = FALSE
      )
    )
    printed <- capture.output(print(result))
    expect_true(any(grepl("Success: 1", printed, fixed = TRUE)))
    expect_true(any(grepl("Results:", printed, fixed = TRUE)))
    expect_true(any(grepl(
      "Resume with: resume_simulation",
      printed,
      fixed = TRUE
    )))
  })
})

describe("F7 end-of-run summary", {
  # cli alerts signal message or warning conditions depending on level;
  # run verbosely, collect every emitted condition into one searchable text
  # block, and return the result alongside it.
  .run_messages <- function(config, ...) {
    conds <- character()
    result <- withCallingHandlers(
      run_simulation(config, ...),
      message = function(m) {
        conds <<- c(conds, conditionMessage(m))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        conds <<- c(conds, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(text = paste(conds, collapse = "\n"), result = result)
  }

  it("prints completion, counts, and path; no resume line when done", {
    path <- withr::local_tempdir()
    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 20L),
      metrics = list(),
      n_replicates = 2L,
      seed = 41L,
      result_path = file.path(path, "complete")
    )
    out <- .run_messages(config, resume = "never", progress = FALSE)

    expect_match(out$text, "Simulation complete: 2/2 tasks succeeded")
    expect_match(out$text, "Results:", fixed = TRUE)
    # Nothing left to resume: no resume instruction in the completed block.
    expect_false(grepl("Resume with", out$text, fixed = TRUE))
    expect_s3_class(out$result, "bayesim_simulation_result")
  })

  it("reports failed counts and per-class detail for a finished run", {
    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "lm"),
      data_generator = function(data_spec, task_ctx) {
        stop("Intentional failure")
      },
      fitter = LinearRegressionFitter(n_draws = 20L),
      metrics = list(),
      n_replicates = 3L,
      seed = 42L,
      max_errors = Inf
    )
    out <- .run_messages(config, resume = "never", progress = FALSE)

    expect_match(out$text, "Simulation finished: 0/3 tasks succeeded, 3 failed")
    expect_match(out$text, "3 task(s) failed", fixed = TRUE)
    # All tasks are terminal: no resume instruction.
    expect_false(grepl("Resume with", out$text, fixed = TRUE))
  })

  it("prints remaining work and the literal resume call after max_errors", {
    path <- withr::local_tempdir()
    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "lm"),
      data_generator = function(data_spec, task_ctx) {
        stop("Intentional failure")
      },
      fitter = LinearRegressionFitter(n_draws = 20L),
      metrics = list(),
      n_replicates = 4L,
      seed = 43L,
      max_errors = 1L,
      result_path = file.path(path, "stopped")
    )
    out <- .run_messages(config, resume = "never", progress = FALSE)

    expect_match(
      out$text,
      "Simulation stopped early \\(error budget exhausted\\): 0 succeeded, 1 failed, 3 of 4 tasks not run"
    )
    # The exact, copyable resume command for this run's path.
    expect_match(out$text, 'Resume with: resume_simulation\\("[^"]+/stopped"')
    expect_false(grepl("Simulation complete", out$text, fixed = TRUE))
    expect_equal(sum(out$result$task_grid$status == "skipped"), 3L)
  })

  it("labels adaptive stops and prints the resume call", {
    path <- withr::local_tempdir()
    .conj_gen <- function(data_spec, task_ctx) {
      n <- data_spec$n
      b <- data_spec$beta
      x <- stats::rnorm(n)
      y <- 1 + b * x + stats::rnorm(n)
      list(
        train = data.frame(y = y, x = x),
        test = NULL,
        response = "y",
        true_params = c(Intercept = 1, x = b, sigma = 1),
        vars_of_interest = c("Intercept", "x", "sigma"),
        meta = list()
      )
    }
    config <- simulation_config(
      data_grid = data.frame(n = 60L, beta = c(0.5, 1)),
      fit_grid = data.frame(model = "lm"),
      data_generator = .conj_gen,
      fitter = LinearRegressionFitter(n_draws = 100L),
      metrics = list(posterior_summary_metric()),
      n_replicates = 6L,
      seed = 44L,
      result_path = file.path(path, "adaptive"),
      checkpoint_every = 2L,
      stop_on = list(
        estimand = "x",
        measure = "bias",
        target_mcse = 100,
        min_reps = 2L,
        check_every = 2L
      )
    )
    out <- .run_messages(config, resume = "never", progress = FALSE)

    expect_match(
      out$text,
      "Simulation stopped early \\(adaptive stopping target reached\\)"
    )
    expect_match(out$text, "tasks not run")
    expect_match(out$text, "Resume with: resume_simulation", fixed = TRUE)
  })

  it("stays silent when verbose = FALSE", {
    path <- withr::local_tempdir()
    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "lm"),
      data_generator = function(data_spec, task_ctx) {
        stop("Intentional failure")
      },
      fitter = LinearRegressionFitter(n_draws = 20L),
      metrics = list(),
      n_replicates = 4L,
      seed = 45L,
      max_errors = 1L,
      result_path = file.path(path, "quiet")
    )
    expect_silent(
      run_simulation(
        config,
        resume = "never",
        progress = FALSE,
        verbose = FALSE
      )
    )
  })

  it("falls back to task results when the grid is absent", {
    tr <- new_task_result(
      task_id = "t1",
      status = "success",
      metrics = list(),
      timing = list(total = 0.1)
    )
    result <- new_simulation_result(
      config_fingerprint = "abc",
      task_results = list(tr),
      timing = list(total = 0.2)
    )
    expect_message(
      print_run_summary(result),
      "Simulation complete: 1/1 tasks succeeded"
    )
  })
})

describe("F6 truthful resume instructions", {
  configless_line <- 'Resume with: resume_simulation\\("[^"]+"\\)$'

  it("keeps the configless command when the manifest is rehydratable", {
    path <- withr::local_tempdir()
    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "lm"),
      data_generator = bayesim:::bayesim_example_data_generator,
      fitter = LinearRegressionFitter(n_draws = 20L),
      metrics = list(),
      n_replicates = 1L,
      seed = 21L,
      result_path = file.path(path, "rehydratable")
    )
    result <- run_simulation(
      config,
      resume = "never",
      progress = FALSE,
      verbose = FALSE
    )
    printed <- capture.output(print(result))
    expect_true(any(grepl(configless_line, printed)))
    expect_false(any(grepl("config =", printed, fixed = TRUE)))
  })

  it("advertises only a configless command that actually works end-to-end", {
    path <- withr::local_tempdir()
    config <- simulation_config(
      data_grid = data.frame(n = 20, beta = 1, sigma = 1),
      fit_grid = data.frame(model = "lm"),
      data_generator = bayesim:::bayesim_example_data_generator,
      fitter = LinearRegressionFitter(n_draws = 20L),
      metrics = list(),
      n_replicates = 1L,
      seed = 31L,
      result_path = file.path(path, "rehydratable-works")
    )
    result <- run_simulation(
      config,
      resume = "never",
      progress = FALSE,
      verbose = FALSE
    )

    status <- run_manifest_rehydration_status(result$checkpoint_path)
    expect_true(status$rehydratable)
    expect_length(status$components, 0L)

    # The advertised command must not merely print; it must succeed.
    resumed <- expect_silent(resume_simulation(
      result$checkpoint_path,
      progress = FALSE,
      verbose = FALSE
    ))
    expect_s3_class(resumed, "bayesim_simulation_result")
  })

  it("asks for config when a closure generator cannot be rehydrated", {
    path <- withr::local_tempdir()
    gen <- fixed_truth_generator(
      truth = c(x = 1),
      draw_data = function(data_spec, task_ctx) {
        x <- stats::rnorm(data_spec$n)
        data.frame(y = x + stats::rnorm(data_spec$n), x = x)
      }
    )
    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "lm"),
      data_generator = gen,
      fitter = LinearRegressionFitter(n_draws = 20L),
      metrics = list(),
      n_replicates = 1L,
      seed = 22L,
      result_path = file.path(path, "closure")
    )
    result <- run_simulation(
      config,
      resume = "never",
      progress = FALSE,
      verbose = FALSE
    )

    printed <- capture.output(print(result))

    # The configless command would fail for this run; it must not be printed.
    expect_false(any(grepl(configless_line, printed)))
    resume_line <- printed[grepl("Resume with", printed)]
    expect_length(resume_line, 1)
    # The printed command must itself be copyable, valid R syntax.
    cmd_expr <- sub("^\\s*Resume with:\\s*", "", resume_line, perl = TRUE)
    parsed <- expect_silent(parse(text = cmd_expr))[[1]]
    expect_identical(parsed[[1]], as.name("resume_simulation"))
    expect_identical(parsed$config, as.name("config"))
    expect_match(deparse1(parsed), "config = config", fixed = TRUE)
    # The reason is actionable: the offending component is named.
    expect_true(any(grepl("data_generator", printed, fixed = TRUE)))
  })

  it("rejects the configless path it no longer advertises for closure runs", {
    path <- withr::local_tempdir()
    gen <- fixed_truth_generator(
      truth = c(x = 1),
      draw_data = function(data_spec, task_ctx) {
        x <- stats::rnorm(data_spec$n)
        data.frame(y = x + stats::rnorm(data_spec$n), x = x)
      }
    )
    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "lm"),
      data_generator = gen,
      fitter = LinearRegressionFitter(n_draws = 20L),
      metrics = list(),
      n_replicates = 1L,
      seed = 32L,
      result_path = file.path(path, "closure-reject")
    )
    result <- run_simulation(
      config,
      resume = "never",
      progress = FALSE,
      verbose = FALSE
    )

    # Pin the failure that motivates the truthful guidance: the exact
    # configless command rehydration used to advertise is rejected.
    expect_error(
      resume_simulation(result$checkpoint_path, progress = FALSE),
      "non-rehydratable"
    )

    # The manifest inspection behind the guidance is component-specific.
    status <- run_manifest_rehydration_status(result$checkpoint_path)
    expect_false(status$rehydratable)
    expect_true("data_generator" %in% status$components)
  })
})
