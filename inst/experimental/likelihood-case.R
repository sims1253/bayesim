# Source-derived design tables for the experimental likelihood vignette.
# Choice-of-likelihood-paper commit 2a96b0d6ac668739363f98f1211c4b0cbc6b9069:
# simulation_study/cont_positive/continuous_positive_simulation.R:5-112
# simulation_study/unit_interval/unit_interval_simulation.R:5-198
# Adaptation: base subset(), explicit link functions, stable names added below.
# These link definitions are stated assumptions, not recovered historical helpers.
likelihood_softplus <- function(x) log(expm1(x))
likelihood_logit <- function(x) stats::qlogis(x)
likelihood_cauchit <- function(x) stats::qcauchy(x)
likelihood_cloglog <- function(x) log(-log1p(-x))

likelihood_positive_design <- function() {
  data_generation_configuration <- expand.grid(
    z1_x_coef = 0.65,
    z1_y_coef = 0.65,
    z2_y_coef = 0.8,
    z3_x_coef = 0.8,
    x_z4_coef = 0.8,
    y_z4_coef = 1,
    sigma_z1 = 0.5,
    sigma_z2 = 0.5,
    sigma_z3 = 0.5,
    sigma_z4 = 0.5,
    sigma_x = 0.5,
    data_N = 100,
    dataset_N = 200,
    data_family = c(
      "gamma",
      "weibull",
      "lognormal",
      "softplusnormal",
      "frechet",
      "inverse.gaussian",
      "betaprime",
      "gompertz"
    ),
    data_link = c("log", "softplus", "identity"),
    lb = 0.000001,
    ub = Inf,
    resample = 1.3,
    x_y_coef = c(NA, 0),
    y_intercept = NA,
    sigma_y = NA,
    shape = c("ramp", "asymmetric", "symmetric"),
    stringsAsFactors = FALSE
  )

  data_generation_configuration <- subset(
    data_generation_configuration,
    !(data_link == "identity" &
      data_family != "lognormal" &
      data_family != "softplusnormal")
  )

  data_generation_configuration <- subset(
    data_generation_configuration,
    !(data_link != "identity" &
      (data_family == "lognormal" | data_family == "softplusnormal"))
  )

  sigma_y_list <- list(
    "gamma" = c(1, 10, 40),
    "weibull" = c(1, 4, 8),
    "lognormal" = c(2, 0.35, 0.15),
    "softplusnormal" = c(1.5, 4, 2),
    "frechet" = c(2, 5, 10),
    "inverse.gaussian" = c(1, 100, 1000),
    "betaprime" = c(1, 10, 50),
    "gompertz" = c(0.2, 0.3, 0.6)
  )

  y_intercept_list <- list(
    "log" = log(c(1, 10, 10)),
    "softplus" = likelihood_softplus(c(1, 10, 10))
  )

  x_y_coef_list <- list(
    "log" = 0.25,
    "softplus" = 0.25
  )

  for (i in seq_len(nrow(data_generation_configuration))) {
    family <- data_generation_configuration$data_family[[i]]
    link <- switch(
      family,
      "lognormal" = "log",
      "softplusnormal" = "softplus",
      data_generation_configuration$data_link[[i]]
    )

    if (is.na(data_generation_configuration$x_y_coef[[i]])) {
      data_generation_configuration$x_y_coef[[i]] <- x_y_coef_list[[link]]
    }
    if (data_generation_configuration$shape[[i]] == "ramp") {
      data_generation_configuration$sigma_y[[i]] <- sigma_y_list[[family]][[1]]
      data_generation_configuration$y_intercept[[i]] <- y_intercept_list[[
        link
      ]][[1]]
    }
    if (data_generation_configuration$shape[[i]] == "asymmetric") {
      data_generation_configuration$sigma_y[[i]] <- sigma_y_list[[family]][[2]]
      data_generation_configuration$y_intercept[[i]] <- y_intercept_list[[
        link
      ]][[2]]
    }
    if (data_generation_configuration$shape[[i]] == "symmetric") {
      data_generation_configuration$sigma_y[[i]] <- sigma_y_list[[family]][[3]]
      data_generation_configuration$y_intercept[[i]] <- y_intercept_list[[
        link
      ]][[3]]
    }
  }
  data_generation_configuration$id <- as.numeric(rownames(
    data_generation_configuration
  ))

  fit_configuration <- expand.grid(
    fit_family = c(
      "inverse.gaussian",
      "gompertz",
      "gamma",
      "weibull",
      "lognormal",
      "softplusnormal",
      "frechet",
      "betaprime",
      "gaussian"
    ),
    fit_link = c("log", "softplus", "identity"),
    formula = c(
      "y ~ x + z1 + z2",
      "y ~ x + z2",
      "y ~ x + z1",
      "y ~ x + z1 + z2 + z3",
      "y ~ x + z1 + z2 + z4"
    ),
    stringsAsFactors = FALSE
  )

  fit_configuration <- subset(
    fit_configuration,
    !(fit_link == "identity" &
      fit_family != "gaussian" &
      fit_family != "lognormal" &
      fit_family != "softplusnormal" &
      fit_family != "lognormal_custom")
  )

  fit_configuration <- subset(
    fit_configuration,
    !(fit_link != "identity" &
      (fit_family == "lognormal" |
        fit_family == "softplusnormal" |
        fit_family == "lognormal_custom"))
  )

  list(
    conditions = data_generation_configuration,
    candidates = fit_configuration
  )
}

likelihood_unit_design <- function() {
  data_generation_configuration <- expand.grid(
    z1_x_coef = -0.45,
    z1_y_coef = 0.45,
    z2_y_coef = 0.6,
    z3_x_coef = 0.8,
    x_z4_coef = 0.45,
    y_z4_coef = 0.8,
    sigma_z1 = 0.5,
    sigma_z2 = 0.5,
    sigma_z3 = 0.5,
    sigma_z4 = 0.5,
    sigma_x = 0.5,
    data_N = 100,
    dataset_N = 200,
    data_family = c(
      "beta",
      "kumaraswamy",
      "logitnormal",
      "cauchitnormal",
      "cloglognormal",
      "simplex"
    ),
    data_link = c("logit", "cauchit", "cloglog", "identity"),
    lb = 0.000001,
    ub = 0.999999,
    resample = 1.3,
    x_y_coef = c(NA, 0),
    y_intercept = NA,
    sigma_y = NA,
    shape = c("symmetric", "asymmetric", "bathtub"),
    stringsAsFactors = FALSE
  )
  data_generation_configuration <- subset(
    data_generation_configuration,
    !(data_link == "identity" &
      data_family != "logitnormal" &
      data_family != "cauchitnormal" &
      data_family != "cloglognormal")
  )

  data_generation_configuration <- subset(
    data_generation_configuration,
    !(data_link != "identity" &
      (data_family == "logitnormal" |
        data_family == "cauchitnormal" |
        data_family == "cloglognormal"))
  )

  sigma_y_list <- list(
    "beta" = c(10, 10, 1.5),
    "kumaraswamy" = c(4, 2.25, 0.5),
    "logitnormal" = c(0.65, 0.8, 2.5),
    "cauchitnormal" = c(0.4, 1, 6),
    "cloglognormal" = c(0.3, 0.4, 3),
    "simplex" = c(1, 1.5, 8)
  )

  y_intercept_list <- list(
    "logit" = c(
      likelihood_logit(0.5),
      likelihood_logit(0.25),
      likelihood_logit(0.5)
    ),
    "cauchit" = c(
      likelihood_cauchit(0.5),
      likelihood_cauchit(0.25),
      likelihood_cauchit(0.5)
    ),
    "cloglog" = c(
      likelihood_cloglog(0.5),
      likelihood_cloglog(0.25),
      likelihood_cloglog(0.5)
    )
  )

  x_y_coef_list <- list(
    "bathtub" = list(
      "logit" = list(
        "beta" = 0.51,
        "kumaraswamy" = 1.22,
        "logitnormal" = 0.98,
        "cauchitnormal" = 0.98,
        "cloglognormal" = 0.98,
        "simplex" = 0.63
      ),
      "cauchit" = list(
        "beta" = 0.59,
        "kumaraswamy" = 1.4,
        "logitnormal" = 2.5,
        "cauchitnormal" = 2.5,
        "cloglognormal" = 2.5,
        "simplex" = 0.58
      ),
      "cloglog" = list(
        "beta" = 0.49,
        "kumaraswamy" = 0.94,
        "logitnormal" = 2,
        "cauchitnormal" = 2,
        "cloglognormal" = 2,
        "simplex" = 0.42
      )
    ),
    "asymmetric" = list(
      "logit" = list(
        "beta" = 0.42,
        "kumaraswamy" = 0.31,
        "logitnormal" = 0.4,
        "cauchitnormal" = 0.4,
        "cloglognormal" = 0.4,
        "simplex" = 0.25
      ),
      "cauchit" = list(
        "beta" = 0.42,
        "kumaraswamy" = 0.44,
        "logitnormal" = 0.5,
        "cauchitnormal" = 0.5,
        "cloglognormal" = 0.5,
        "simplex" = 0.31
      ),
      "cloglog" = list(
        "beta" = 0.35,
        "kumaraswamy" = 0.3,
        "logitnormal" = 0.2,
        "cauchitnormal" = 0.2,
        "cloglognormal" = 0.2,
        "simplex" = 0.21
      )
    ),
    "symmetric" = list(
      "logit" = list(
        "beta" = 0.32,
        "kumaraswamy" = 0.26,
        "logitnormal" = 0.4,
        "cauchitnormal" = 0.4,
        "cloglognormal" = 0.4,
        "simplex" = 0.3
      ),
      "cauchit" = list(
        "beta" = 0.28,
        "kumaraswamy" = 0.25,
        "logitnormal" = 0.23,
        "cauchitnormal" = 0.23,
        "cloglognormal" = 0.23,
        "simplex" = 0.26
      ),
      "cloglog" = list(
        "beta" = 0.25,
        "kumaraswamy" = 0.20,
        "logitnormal" = 0.18,
        "cauchitnormal" = 0.18,
        "cloglognormal" = 0.18,
        "simplex" = 0.17
      )
    )
  )

  bathtub_list <- list(
    z1_x_coef = -0.85,
    z1_y_coef = 0.85,
    z2_y_coef = 0.6,
    z3_x_coef = 0.8,
    x_z4_coef = 0.35,
    y_z4_coef = 0.5
  )

  for (i in seq_len(nrow(data_generation_configuration))) {
    family <- data_generation_configuration$data_family[[i]]
    shape <- data_generation_configuration$shape[[i]]
    link <- data_generation_configuration$data_link[[i]]
    link <- switch(
      family,
      "logitnormal" = "logit",
      "cauchitnormal" = "cauchit",
      "cloglognormal" = "cloglog",
      data_generation_configuration$data_link[[i]]
    )

    if (is.na(data_generation_configuration$x_y_coef[[i]])) {
      data_generation_configuration$x_y_coef[[i]] <- x_y_coef_list[[shape]][[
        link
      ]][[family]]
    }
    if (shape == "symmetric") {
      data_generation_configuration$sigma_y[[i]] <- sigma_y_list[[family]][[1]]
      data_generation_configuration$y_intercept[[i]] <- y_intercept_list[[
        link
      ]][[1]]
    }
    if (shape == "asymmetric") {
      data_generation_configuration$sigma_y[[i]] <- sigma_y_list[[family]][[2]]
      data_generation_configuration$y_intercept[[i]] <- y_intercept_list[[
        link
      ]][[2]]
    }
    if (shape == "bathtub") {
      data_generation_configuration$sigma_y[[i]] <- sigma_y_list[[family]][[3]]
      data_generation_configuration$y_intercept[[i]] <- y_intercept_list[[
        link
      ]][[3]]
      data_generation_configuration$z1_x_coef[[i]] <- bathtub_list$z1_x_coef
      data_generation_configuration$z1_y_coef[[i]] <- bathtub_list$z1_y_coef
      data_generation_configuration$z2_y_coef[[i]] <- bathtub_list$z2_y_coef
      data_generation_configuration$z3_x_coef[[i]] <- bathtub_list$z3_x_coef
      data_generation_configuration$x_z4_coef[[i]] <- bathtub_list$x_z4_coef
      data_generation_configuration$y_z4_coef[[i]] <- bathtub_list$y_z4_coef
    }
  }
  data_generation_configuration$id <- as.numeric(rownames(
    data_generation_configuration
  ))

  fit_configuration <- expand.grid(
    fit_family = c(
      "beta",
      "kumaraswamy",
      "logitnormal",
      "cauchitnormal",
      "cloglognormal",
      "simplex",
      "gaussian"
    ),
    fit_link = c("logit", "cauchit", "cloglog", "identity"),
    formula = c(
      "y ~ x + z1 + z2",
      "y ~ x + z2",
      "y ~ x + z1",
      "y ~ x + z1 + z2 + z3",
      "y ~ x + z1 + z2 + z4"
    ),
    stringsAsFactors = FALSE
  )

  fit_configuration <- subset(
    fit_configuration,
    !(fit_link == "identity" &
      fit_family != "gaussian" &
      fit_family != "logitnormal" &
      fit_family != "cauchitnormal" &
      fit_family != "cloglognormal")
  )

  fit_configuration <- subset(
    fit_configuration,
    !(fit_link != "identity" &
      (fit_family == "logitnormal" |
        fit_family == "cauchitnormal" |
        fit_family == "cloglognormal"))
  )

  list(
    conditions = data_generation_configuration,
    candidates = fit_configuration
  )
}

likelihood_effective_link <- function(family, link) {
  replacements <- c(
    lognormal = "log",
    softplusnormal = "softplus",
    logitnormal = "logit",
    cauchitnormal = "cauchit",
    cloglognormal = "cloglog"
  )
  found <- family %in% names(replacements)
  link[found] <- unname(replacements[family[found]])
  link
}

likelihood_design <- function(domain = c("positive", "unit")) {
  domain <- match.arg(domain)
  design <- if (domain == "positive") {
    likelihood_positive_design()
  } else {
    likelihood_unit_design()
  }
  d <- design$conditions
  d$condition_id <- paste(
    domain,
    d$data_family,
    d$data_link,
    d$shape,
    ifelse(d$x_y_coef == 0, "null", "effect"),
    sep = ":"
  )
  d$effective_link <- likelihood_effective_link(d$data_family, d$data_link)
  m <- design$candidates
  m$effective_link <- likelihood_effective_link(m$fit_family, m$fit_link)
  m$method_id <- paste(m$fit_family, m$fit_link, m$formula, sep = ":")
  stopifnot(
    !anyDuplicated(d$condition_id),
    !anyDuplicated(m$method_id),
    !anyNA(d$x_y_coef),
    !anyNA(d$y_intercept),
    !anyNA(d$sigma_y)
  )
  list(conditions = d, candidates = m)
}

# Illustrative SEM with source coefficients. This is not the historical
# generator: it omits response-bound rejection/resampling and supports only
# Gamma/log and Beta/logit. Unsupported conditions fail instead of changing DGP.
likelihood_generate_subset <- function(condition, context) {
  c <- condition
  supported <- (c$data_family == "gamma" && c$data_link == "log") ||
    (c$data_family == "beta" && c$data_link == "logit")
  if (!supported) {
    stop(
      "The illustrative generator supports only Gamma/log and Beta/logit; supply the historical generator for other conditions."
    )
  }
  draw <- function() {
    n <- c$data_N
    z1 <- stats::rnorm(n, sd = c$sigma_z1)
    z2 <- stats::rnorm(n, sd = c$sigma_z2)
    z3 <- stats::rnorm(n, sd = c$sigma_z3)
    x <- c$z1_x_coef * z1 + c$z3_x_coef * z3 + stats::rnorm(n, sd = c$sigma_x)
    eta <- c$y_intercept + c$x_y_coef * x + c$z1_y_coef * z1 + c$z2_y_coef * z2
    y <- if (c$data_family == "gamma") {
      mu <- exp(eta)
      stats::rgamma(n, shape = c$sigma_y, rate = c$sigma_y / mu)
    } else {
      mu <- stats::plogis(eta)
      stats::rbeta(n, shape1 = mu * c$sigma_y, shape2 = (1 - mu) * c$sigma_y)
    }
    z4 <- c$x_z4_coef * x + c$y_z4_coef * y + stats::rnorm(n, sd = c$sigma_z4)
    data.frame(y, x, z1, z2, z3, z4)
  }
  list(
    train = draw(),
    test = draw(),
    truth = c$x_y_coef,
    effective_link = c$effective_link
  )
}

likelihood_fit <- function(data, context) {
  settings <- context$settings
  # Built-in brms families only. The full source grid also contains custom
  # families. brms reports unsupported configurations when execution is tried.
  brms::brm(
    formula = stats::as.formula(settings$formula),
    data = data$train,
    family = brms::brmsfamily(settings$fit_family, link = settings$fit_link),
    backend = "cmdstanr",
    seed = context$seed,
    chains = 4L,
    cores = 1L,
    iter = 2000L,
    refresh = 0
  )
}

likelihood_extract_draws <- function(fit, data, context) {
  posterior::subset_draws(posterior::as_draws_array(fit), variable = "b_x")
}

likelihood_extract_diagnostics <- function(fit, data, context) {
  summary <- posterior::summarise_draws(posterior::as_draws_array(fit))
  nuts <- brms::nuts_params(fit)
  list(
    rhat = max(summary$rhat, na.rm = TRUE),
    ess_bulk = min(summary$ess_bulk, na.rm = TRUE),
    ess_tail = min(summary$ess_tail, na.rm = TRUE),
    divergents = sum(nuts$Value[nuts$Parameter == "divergent__"])
  )
}

likelihood_extract_pointwise <- function(fit, data, context) {
  brms::loo(fit)$pointwise
}

likelihood_measure <- function(artifacts, context) {
  draws <- as.numeric(artifacts$draws)
  truth <- artifacts$data$truth
  # This mapping is an explicit convention for this illustrative SEM. Other
  # studies must supply their own estimand mapping, even if names/links match.
  comparable <- context$settings$effective_link == artifacts$data$effective_link
  interval <- stats::quantile(draws, c(0.025, 0.975))
  diagnostics <- artifacts$diagnostics
  data.frame(
    target = "x",
    value = if (comparable) mean(draws) - truth else NA_real_,
    comparable = comparable,
    rmse_s = if (comparable) sqrt(mean((draws - truth)^2)) else NA_real_,
    sign_probability = mean(draws > 0),
    sig95 = interval[1] > 0 || interval[2] < 0,
    elpd_loo = sum(artifacts$pointwise[, "elpd_loo"]),
    rhat = diagnostics$rhat,
    ess_bulk = diagnostics$ess_bulk,
    ess_tail = diagnostics$ess_tail,
    divergents = diagnostics$divergents,
    row.names = NULL
  )
}

# One group is a dataset and a candidate formula. Raw prediction comparisons
# include every finite score, regardless of coefficient comparability.
likelihood_elpd_gap <- function(results, artifacts, context) {
  x <- results[results$measurement == "performance", , drop = FALSE]
  valid <- is.finite(x$elpd_loo)
  if (!any(valid)) {
    return(data.frame(
      candidate = character(),
      value = numeric(),
      n_reference = integer()
    ))
  }
  x <- x[valid, , drop = FALSE]
  data.frame(
    candidate = x$method_id,
    value = x$elpd_loo - max(x$elpd_loo),
    n_reference = nrow(x)
  )
}

# Explicit recovery-analysis policy, inspired by the source filters. The
# matching-effective-link restriction is specific to this estimand mapping.
likelihood_selection <- function(measurements, attempts) {
  x <- measurements
  eligible <- x$measurement == "performance" &
    x$comparable &
    is.finite(x$elpd_loo) &
    is.finite(x$value) &
    is.finite(x$rmse_s) &
    is.finite(x$rhat) &
    is.finite(x$ess_bulk) &
    is.finite(x$ess_tail) &
    is.finite(x$divergents) &
    x$divergents <= 10 &
    x$rhat <= 1.01 &
    x$ess_bulk > 400 &
    x$ess_tail > 400 &
    abs(x$value) < 10 &
    x$rmse_s < 10
  x[!is.na(eligible) & eligible, , drop = FALSE]
}

likelihood_declaration <- function(design, name) {
  s <- study(
    name,
    conditions = design$conditions,
    generate = likelihood_generate_subset,
    version = "illustrative-sem-1"
  )
  methods <- lapply(seq_len(nrow(design$candidates)), function(i) {
    settings <- as.list(design$candidates[i, , drop = FALSE])
    study_method(
      likelihood_fit,
      extract = list(
        draws = likelihood_extract_draws,
        diagnostics = likelihood_extract_diagnostics,
        pointwise = likelihood_extract_pointwise
      ),
      settings = settings,
      version = "brms-subset-1"
    )
  })
  names(methods) <- design$candidates$method_id
  s <- do.call(with_methods, c(list(s), methods))
  s <- with_measure(
    s,
    "performance",
    likelihood_measure,
    needs = c("data", "draws", "diagnostics", "pointwise")
  )
  s <- with_comparison(
    s,
    "elpd_gap",
    likelihood_elpd_gap,
    by = c("dataset_id", "method_formula")
  )
  with_retention(s, artifacts = c("draws", "diagnostics", "pointwise"))
}
