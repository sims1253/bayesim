# I3: adaptive stopping on MCSE targets.
library(bayesim)

# Conjugate linear regression generator (mirrors test-performance-measures.R).
.gen <- function(data_spec, task_ctx) {
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

describe("adaptive stopping (stop_on)", {
  it("requires every condition cell to meet the MCSE target", {
    task_results <- list(new_task_result(
      task_id = "t1",
      status = "success",
      metrics = list(),
      timing = list(total = 0)
    ))
    task_grid <- data.frame(
      task_id = c("t1", "t2"),
      data_idx = c(1L, 2L),
      fit_idx = c(1L, 1L),
      status = c("success", "success")
    )
    stop_on <- list(
      estimand = "x",
      measure = "bias",
      target_mcse = 0.1,
      min_reps = 2L
    )

    local_mocked_bindings(
      bayesim_adaptive_summary = function(...) data.frame(dummy = 1),
      performance_measures = function(...) {
        data.frame(
          condition = c("stable", "noisy"),
          measure = c("bias", "bias"),
          mcse = c(0.01, 0.2),
          n_sim = c(2L, 2L)
        )
      },
      .package = "bayesim"
    )

    expect_false(bayesim_adaptive_check(
      task_results,
      task_grid,
      stop_on,
      config = NULL
    ))
  })

  it("stops only when every selected MCSE is finite and below target", {
    task_results <- list(new_task_result(
      task_id = "t1",
      status = "success",
      metrics = list(),
      timing = list(total = 0)
    ))
    task_grid <- data.frame(
      task_id = "t1",
      data_idx = 1L,
      fit_idx = 1L,
      status = "success"
    )
    stop_on <- list(
      estimand = "x",
      measure = "bias",
      target_mcse = 0.1,
      min_reps = 2L
    )

    local_mocked_bindings(
      bayesim_adaptive_summary = function(...) data.frame(dummy = 1),
      performance_measures = function(...) {
        data.frame(measure = "bias", mcse = 0.01, n_sim = 2L)
      },
      .package = "bayesim"
    )

    expect_true(bayesim_adaptive_check(
      task_results,
      task_grid,
      stop_on,
      config = NULL
    ))
  })

  it("stops early when a loose target_mcse is met", {
    # Loose target_mcse = 0.5 on bias should be satisfied after a handful of
    # reps. min_reps=2, check_every=2 => the check first fires at n_completed=2.
    cfg <- simulation_config(
      data_grid = data.frame(n = 60L, beta = 0.5),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 100L),
      metrics = list(posterior_summary_metric()),
      n_replicates = 20L,
      seed = 11L,
      checkpoint_every = 2L,
      stop_on = list(
        estimand = "x",
        measure = "bias",
        target_mcse = 0.5,
        min_reps = 2L,
        check_every = 2L
      )
    )
    res <- run_simulation(cfg, resume = "never", progress = FALSE)
    summ <- res$summary

    n_success <- sum(summ$status == "success", na.rm = TRUE)
    n_skipped <- sum(summ$status == "skipped", na.rm = TRUE)
    # Either the run stopped short (skipped tasks present) or, if the MCSE
    # check happened to never trip, all 20 ran. The loose target should trip.
    expect_true(n_success < 20L || n_skipped > 0L)
    # At least min_reps tasks ran.
    expect_gte(n_success, 2L)
    # Skipped tasks appear in the summary (status machinery supports it).
    if (n_success < 20L) {
      expect_gt(n_skipped, 0L)
    }
  })

  it("does not stop before every condition has min_reps usable estimates", {
    cfg <- simulation_config(
      data_grid = data.frame(n = 60L, beta = c(0.5, 1)),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 100L),
      metrics = list(posterior_summary_metric()),
      n_replicates = 4L,
      seed = 12L,
      checkpoint_every = 2L,
      stop_on = list(
        estimand = "x",
        measure = "bias",
        target_mcse = 100,
        min_reps = 2L,
        check_every = 2L
      )
    )

    res <- run_simulation(cfg, resume = "never", progress = FALSE)
    successes <- table(factor(
      res$summary$data_beta[res$summary$status == "success"],
      levels = c(0.5, 1)
    ))

    expect_true(all(successes >= 2L))
    expect_gt(sum(res$summary$status == "skipped"), 0L)
  })

  it("persists the precision snapshot and replicate round that triggered stop", {
    result_path <- file.path(withr::local_tempdir(), "adaptive-state")
    cfg <- simulation_config(
      data_grid = data.frame(n = 60L, beta = c(0.5, 1)),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 100L),
      metrics = list(posterior_summary_metric()),
      n_replicates = 4L,
      seed = 18L,
      result_path = result_path,
      checkpoint_every = 2L,
      stop_on = list(
        estimand = "x",
        measure = "bias",
        target_mcse = 100,
        min_reps = 2L,
        check_every = 2L
      )
    )

    run_simulation(cfg, resume = "never", progress = FALSE, verbose = FALSE)
    checkpoint <- read_checkpoint(result_path)
    state <- checkpoint$adaptive_state

    expect_true(isTRUE(state$triggered))
    expect_equal(state$estimand, "x")
    expect_equal(state$measure, "bias")
    expect_gte(state$completed_rounds, 2L)
    expect_equal(length(state$cells), 2L)
    expect_true(all(vapply(
      state$cells,
      function(x) is.finite(x$mcse),
      logical(1)
    )))
  })

  it("checks at replicate rounds independently of checkpoint cadence", {
    cfg <- simulation_config(
      data_grid = data.frame(n = 60L, beta = c(0.5, 1)),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 100L),
      metrics = list(posterior_summary_metric()),
      n_replicates = 30L,
      seed = 28L,
      checkpoint_every = 50L,
      stop_on = list(
        estimand = "x",
        measure = "bias",
        target_mcse = 100,
        min_reps = 2L,
        check_every = 2L
      )
    )

    result <- run_simulation(
      cfg,
      resume = "never",
      progress = FALSE,
      verbose = FALSE
    )
    successes <- table(factor(
      result$summary$data_beta[result$summary$status == "success"],
      levels = c(0.5, 1)
    ))

    expect_equal(as.integer(successes), c(2L, 2L))
  })

  it("NULL stop_on runs every task (control)", {
    cfg <- simulation_config(
      data_grid = data.frame(n = 60L, beta = 0.5),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 100L),
      metrics = list(posterior_summary_metric()),
      n_replicates = 10L,
      seed = 23L
      # stop_on omitted => NULL
    )
    res <- run_simulation(cfg, resume = "never", progress = FALSE)
    summ <- res$summary
    expect_equal(sum(summ$status == "success", na.rm = TRUE), 10L)
    expect_equal(sum(summ$status == "skipped", na.rm = TRUE), 0L)
  })

  it("rejects a malformed stop_on policy", {
    expect_error(
      simulation_config(
        data_grid = data.frame(n = 60L, beta = 0.5),
        fit_grid = data.frame(model = "lm"),
        data_generator = .gen,
        fitter = LinearRegressionFitter(n_draws = 10L),
        metrics = list(posterior_summary_metric()),
        n_replicates = 2L,
        seed = 1L,
        stop_on = list(estimand = "x", measure = "bias") # missing target_mcse
      ),
      class = "bayesim_config_error"
    )
    expect_error(
      simulation_config(
        data_grid = data.frame(n = 60L, beta = 0.5),
        fit_grid = data.frame(model = "lm"),
        data_generator = .gen,
        fitter = LinearRegressionFitter(n_draws = 10L),
        metrics = list(posterior_summary_metric()),
        n_replicates = 2L,
        seed = 1L,
        stop_on = list(estimand = "x", measure = "bogus", target_mcse = 0.1)
      ),
      class = "bayesim_config_error"
    )
  })
})
