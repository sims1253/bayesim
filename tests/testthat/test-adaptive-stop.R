# I3: adaptive stopping on MCSE targets.
skip_on_cran()
library(bayesim)

# Conjugate linear regression generator (mirrors test-performance-measures.R).
.gen <- function(data_spec, task_ctx) {
  n <- data_spec$n; b <- data_spec$beta
  x <- stats::rnorm(n)
  y <- 1 + b * x + stats::rnorm(n)
  list(
    train = data.frame(y = y, x = x), test = NULL, response = "y",
    true_params = c(Intercept = 1, x = b, sigma = 1),
    vars_of_interest = c("Intercept", "x", "sigma"),
    meta = list()
  )
}

describe("adaptive stopping (stop_on)", {
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
        estimand = "x", measure = "bias",
        target_mcse = 0.5, min_reps = 2L, check_every = 2L
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
        n_replicates = 2L, seed = 1L,
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
        n_replicates = 2L, seed = 1L,
        stop_on = list(estimand = "x", measure = "bogus", target_mcse = 0.1)
      ),
      class = "bayesim_config_error"
    )
  })
})
