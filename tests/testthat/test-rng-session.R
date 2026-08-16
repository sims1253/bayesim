library(bayesim)

.rng_gen <- function(data_spec, task_ctx) {
  x <- stats::rnorm(data_spec$n)
  list(
    train = data.frame(y = x, x = x),
    test = NULL,
    response = "y",
    true_params = c(Intercept = 0, x = 1, sigma = 0.1),
    vars_of_interest = c("Intercept", "x", "sigma")
  )
}

describe("run_simulation RNG isolation", {
  it("restores the caller's RNG kind and state", {
    original_kind <- RNGkind()
    original_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit(
      {
        do.call(base::RNGkind, as.list(original_kind))
        if (is.null(original_seed)) {
          if (exists(".Random.seed", envir = .GlobalEnv)) {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        } else {
          assign(".Random.seed", original_seed, envir = .GlobalEnv)
        }
      },
      add = TRUE
    )

    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(9182)
    before_kind <- RNGkind()
    before_seed <- get(".Random.seed", envir = .GlobalEnv)

    config <- simulation_config(
      data_grid = data.frame(n = 10L),
      fit_grid = data.frame(model = "linear"),
      data_generator = .rng_gen,
      fitter = LinearRegressionFitter(n_draws = 10L),
      metrics = list(),
      n_replicates = 1L,
      seed = 42L
    )
    run_simulation(config, resume = "never", progress = FALSE)

    expect_identical(RNGkind(), before_kind)
    expect_identical(get(".Random.seed", envir = .GlobalEnv), before_seed)
  })
})
