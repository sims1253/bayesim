# Acceptance tests for the purrr + mirai transport.
# - run_task_safe is total; fatal conditions raised inside a task under
#   daemons stop the run with the original condition class.
# - determinism (sequential == daemons(2)) still holds on the new transport.
# - workers = 2 matches the sequential summary and leaves daemons unset.

.gen <- function(data_spec, task_ctx) {
  n <- data_spec$n %||% 20L
  list(
    train = data.frame(y = stats::rnorm(n), x = stats::rnorm(n)),
    test = NULL,
    response = "y",
    true_params = c(beta = 0),
    vars_of_interest = "beta",
    meta = list()
  )
}

describe("purrr/mirai transport", {
  it("releases model banks on reused daemons after success and fatal failure", {
    mirai::daemons(1)
    on.exit(mirai::daemons(0), add = TRUE)
    withr::local_options(bayesim.model_bank = list(sentinel = "current"))

    gen <- function(data_spec, task_ctx) {
      if (
        !identical(getOption("bayesim.model_bank"), list(sentinel = "current"))
      ) {
        stop(bayesim::bayesim_config_error("model bank was not installed"))
      }
      if (data_spec$fail) {
        stop(bayesim::bayesim_config_error(
          "deliberate failure with model bank"
        ))
      }
      list(train = data.frame(y = 1:5, x = 1:5), response = "y")
    }

    for (fail in c(FALSE, TRUE)) {
      config <- simulation_config(
        data_grid = data.frame(fail = fail),
        fit_grid = data.frame(model = "baseline"),
        data_generator = gen,
        fitter = MockFitter(),
        metrics = list(),
        n_replicates = 3L,
        checkpoint_every = 1L,
        seed = 42L
      )
      config_spec <- as_config_spec(config)
      config_spec$data_generator <- gen
      config_spec$package_name <- "bayesim"
      run <- function() {
        execute_tasks(
          task_grid = create_task_grid(config),
          config = config,
          config_spec = config_spec,
          fitter = config@fitter,
          metrics = config@metrics,
          retain = "metrics",
          max_errors = Inf,
          progress = FALSE,
          verbose = FALSE,
          checkpoint_every = 1L
        )
      }
      if (fail) {
        expect_error(
          run(),
          "deliberate failure with model bank",
          class = "bayesim_config_error"
        )
      } else {
        expect_identical(run()$task_grid$status, rep("success", 3L))
      }
      expect_true(mirai::daemons_set())
      bank <- mirai::call_mirai(mirai::mirai(getOption(
        "bayesim.model_bank"
      )))$data
      expect_null(bank)
    }
  })

  it("clears a stale daemon bank before a study without precompiled models", {
    mirai::daemons(1)
    on.exit(mirai::daemons(0), add = TRUE)
    mirai::everywhere(options(bayesim.model_bank = list(stale = TRUE)))
    config <- simulation_config(
      data_grid = data.frame(n = 5L),
      fit_grid = data.frame(model = "baseline"),
      data_generator = function(data_spec, task_ctx) {
        if (!is.null(getOption("bayesim.model_bank"))) {
          stop(bayesim::bayesim_config_error("stale bank reached the task"))
        }
        list(train = data.frame(y = 1:5, x = 1:5), response = "y")
      },
      fitter = MockFitter(),
      metrics = list(),
      n_replicates = 1L,
      seed = 42L
    )
    result <- run_simulation(config, progress = FALSE, verbose = FALSE)
    expect_identical(result$summary$status, "success")
  })
  it("fatal conditions raised inside a task stop the run under daemons", {
    # A data generator that raises a fatal bayesim_config_error. Generators are
    # crated into the task transport (config_spec$data_generator), so any helper
    # they call must be namespace-qualified to resolve on daemons (bayesim is
    # installed there; bayesim_config_error is exported).
    fatal_gen <- function(data_spec, task_ctx) {
      stop(bayesim::bayesim_config_error(
        "deliberate fatal failure inside a task"
      ))
    }

    config <- simulation_config(
      data_grid = data.frame(n = 20),
      fit_grid = data.frame(model = "baseline"),
      data_generator = fatal_gen,
      fitter = MockFitter(),
      metrics = list(),
      n_replicates = 2L,
      seed = 42L
    )

    # Sanity: sequential also raises.
    expect_error(
      run_simulation(config, resume = "never", progress = FALSE),
      class = "bayesim_config_error"
    )

    mirai::daemons(2)
    on.exit(mirai::daemons(0), add = TRUE)
    expect_error(
      run_simulation(config, resume = "never", progress = FALSE),
      class = "bayesim_config_error"
    )
  })

  it("sequential == daemons(2) summaries match", {
    config <- simulation_config(
      data_grid = data.frame(n = c(30, 60)),
      fit_grid = data.frame(model = "baseline"),
      data_generator = .gen,
      fitter = MockFitter(),
      metrics = list(pred_rmse_metric()),
      n_replicates = 2L,
      seed = 42L
    )
    seq_res <- run_simulation(config, resume = "never", progress = FALSE)
    mirai::daemons(2)
    on.exit(mirai::daemons(0), add = TRUE)
    par_res <- run_simulation(config, resume = "never", progress = FALSE)

    norm <- function(x) {
      x <- x[order(x$task_id), , drop = FALSE]
      x$timing_total <- NULL
      x
    }
    expect_equal(norm(seq_res$summary), norm(par_res$summary))
  })
})

describe("workers convenience argument", {
  it("workers = 2 matches the sequential summary and tears down daemons", {
    config <- simulation_config(
      data_grid = data.frame(n = 30),
      fit_grid = data.frame(model = "baseline"),
      data_generator = .gen,
      fitter = MockFitter(),
      metrics = list(pred_rmse_metric()),
      n_replicates = 2L,
      seed = 42L
    )
    seq_res <- run_simulation(config, resume = "never", progress = FALSE)

    expect_false(mirai::daemons_set())
    par_res <- run_simulation(
      config,
      resume = "never",
      progress = FALSE,
      workers = 2
    )
    expect_false(mirai::daemons_set())

    norm <- function(x) {
      x <- x[order(x$task_id), , drop = FALSE]
      x$timing_total <- NULL
      x
    }
    expect_equal(norm(seq_res$summary), norm(par_res$summary))
  })

  it("errors when workers is non-NULL and daemons are already set", {
    mirai::daemons(2)
    on.exit(mirai::daemons(0), add = TRUE)
    config <- simulation_config(
      data_grid = data.frame(n = 30),
      fit_grid = data.frame(model = "baseline"),
      data_generator = .gen,
      fitter = MockFitter(),
      metrics = list(),
      n_replicates = 1L,
      seed = 42L
    )
    expect_error(
      run_simulation(config, resume = "never", progress = FALSE, workers = 2),
      class = "bayesim_config_error"
    )
  })
})
