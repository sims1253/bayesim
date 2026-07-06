# C1/C2 acceptance tests for the purrr + mirai transport.
# - C1: run_task_safe is total; fatal conditions raised inside a task under
#   daemons stop the run with the original condition class.
# - C1: determinism (sequential == daemons(2)) still holds on the new transport.
# - C2: workers = 2 matches the sequential summary and leaves daemons unset.

skip_on_cran()

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

describe("C1: purrr/mirai transport", {
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

describe("C2: workers convenience argument", {
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
