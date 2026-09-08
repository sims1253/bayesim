# The experimental API is sourced explicitly, just as in its vignettes.
grammar <- new.env(parent = globalenv())
sys.source(
  system.file("experimental", "study-grammar.R", package = "bayesim"),
  grammar
)

make_grammar_study <- function() {
  s <- grammar$study(
    "test",
    data.frame(condition_id = c("small", "large"), n = c(8L, 12L)),
    function(condition, context) {
      stopifnot(is.null(context$method_id), is.null(context$settings))
      rnorm(condition$n)
    }
  )
  method <- grammar$study_method(
    function(data, context) {
      list(
        data_mean = mean(data),
        draws = array(
          rnorm(40, mean(data)),
          c(20, 2, 1),
          dimnames = list(NULL, NULL, "theta")
        )
      )
    },
    extract = list(
      draws = function(fit, data, context) fit$draws,
      truth = function(fit, data, context) fit$data_mean
    )
  )
  s <- grammar$with_methods(s, first = method, second = method)
  s <- grammar$with_measure(
    s,
    "estimate",
    function(artifacts, context) {
      data.frame(
        target = "theta",
        value = mean(artifacts$draws),
        generated = artifacts$truth
      )
    },
    needs = c("draws", "truth")
  )
  grammar$with_retention(s, c("draws", "truth"))
}

canonical_grammar_rows <- function(rows) {
  rows <- rows[
    order(rows$dataset_id, rows$method_id, rows$measurement),
    ,
    drop = FALSE
  ]
  rownames(rows) <- NULL
  rows
}

test_that("plans validate artifact requirements without executing scientific code", {
  s <- grammar$study(
    "plan",
    data.frame(condition_id = "a"),
    function(condition, context) stop("generated"),
    prepare = function(condition, context) stop("prepared")
  )
  s <- grammar$with_methods(
    s,
    method = grammar$study_method(function(data, context) stop("fitted"))
  )
  p <- grammar$plan_study(s, 2000000)
  expect_equal(p$fits, 2000000)
  expect_true(p$prepare)
  missing <- grammar$with_measure(
    s,
    "interval",
    function(artifacts, context) data.frame(),
    needs = "draws"
  )
  expect_error(grammar$plan_study(missing, 1), "cannot supply: draws")
  expect_error(
    grammar$with_comparison(
      s,
      "across",
      function(results, artifacts, context) data.frame(),
      by = "condition_id",
      needs = "fit"
    ),
    "within dataset_id"
  )
  expect_error(grammar$plan_study(s, 0), "replicates")
})

test_that("dataset identities survive method and replicate extensions and condition reordering", {
  s <- make_grammar_study()
  first <- grammar$run_study(s, 2, seed = 102L)
  expect_equal(
    first$measurements$generated[c(1, 3, 5, 7)],
    first$measurements$generated[c(2, 4, 6, 8)]
  )
  expect_false(identical(
    first$measurements$value[c(1, 3, 5, 7)],
    first$measurements$value[c(2, 4, 6, 8)]
  ))
  extended <- grammar$with_methods(s, third = s$methods$first)
  extended$conditions <- extended$conditions[2:1, , drop = FALSE]
  extended$methods <- extended$methods[c("third", "second", "first")]
  later <- grammar$run_study(extended, 3, seed = 102L)
  old <- later$measurements[
    later$measurements$replicate <= 2 & later$measurements$method_id != "third",
  ]
  expect_equal(
    canonical_grammar_rows(first$measurements),
    canonical_grammar_rows(old)
  )
})

test_that("saved draws support additional measurements without refitting", {
  s <- make_grammar_study()
  path <- file.path(withr::local_tempdir(), "run")
  first <- grammar$run_study(s, 2, path = path)
  s <- grammar$with_measure(
    s,
    "interval90",
    function(artifacts, context) {
      stopifnot(identical(dim(artifacts$draws), c(20L, 2L, 1L)))
      q <- quantile(artifacts$draws, c(0.05, 0.95))
      data.frame(target = "theta", lower = q[[1]], upper = q[[2]])
    },
    needs = "draws"
  )
  later <- grammar$run_study(s, 3, path = path)
  old <- later$attempts[later$attempts$replicate <= 2, ]
  expect_equal(first$attempts, old) # persisted timings would change on refit
  expect_equal(sum(later$measurements$measurement == "interval90"), 12)
  expect_equal(
    first$measurements,
    later$measurements[
      later$measurements$replicate <= 2 &
        later$measurements$measurement == "estimate",
      names(first$measurements),
      drop = FALSE
    ],
    ignore_attr = TRUE
  )
  expect_error(grammar$run_study(s, 2, seed = 2, path = path), "seed")
  original <- make_grammar_study()
  no_draws <- grammar$with_retention(original)
  dropped_path <- file.path(withr::local_tempdir(), "discarded")
  grammar$run_study(no_draws, 1, path = dropped_path)
  s$retention <- character()
  expect_error(
    grammar$run_study(s, 1, path = dropped_path),
    "discarded artifacts: draws"
  )
})

test_that("preparation is shared and independently addressable from an external scheduler", {
  count_file <- tempfile()
  s <- make_grammar_study()
  s$conditions$counter <- count_file
  s$prepare <- function(condition, context) {
    cat(condition$condition_id, "\n", file = condition$counter, append = TRUE)
    rnorm(2)
  }
  path <- file.path(withr::local_tempdir(), "prepared")
  first <- grammar$run_study(s, 2, path = path)
  expect_length(readLines(count_file), 2)
  one <- grammar$evaluate_replicate(s, "small", 1, path = path)
  expect_length(readLines(count_file), 2)
  expect_equal(
    one$measurements,
    first$measurements[
      first$measurements$condition_id == "small" &
        first$measurements$replicate == 1,
    ]
  )
  unlink(count_file)
})

test_that("group inputs survive until comparison and complete comparisons are reusable", {
  s <- make_grammar_study()
  s <- grammar$with_comparison(
    s,
    "paired",
    function(results, artifacts, context) {
      stopifnot(length(artifacts) == 2)
      data.frame(value = mean(artifacts$first$draws - artifacts$second$draws))
    },
    needs = "draws"
  )
  s <- grammar$with_retention(s)
  path <- file.path(withr::local_tempdir(), "group")
  first <- grammar$run_study(s, 2, path = path)
  expect_equal(nrow(first$comparisons), 4)
  second <- grammar$run_study(s, 2, path = path)
  expect_equal(first$comparisons, second$comparisons)
  expect_equal(first$attempts, second$attempts)
  s$methods <- rev(s$methods)
  reordered <- grammar$run_study(s, 2, path = path)
  expect_equal(first$comparisons, reordered$comparisons)
})

test_that("analysis policies update reference groups without hiding failed methods", {
  s <- make_grammar_study()
  s <- grammar$with_methods(
    s,
    failed = grammar$study_method(
      function(data, context) stop("deliberate failure"),
      extract = s$methods$first$extract
    )
  )
  s <- grammar$with_comparison(s, "gap", function(results, artifacts, context) {
    good <- results[results$measurement == "estimate", , drop = FALSE]
    data.frame(
      n_attempted = nrow(context$attempts),
      n_failed = sum(context$attempts$status == "failed"),
      n_used = nrow(good),
      best = if (nrow(good)) max(good$value) else NA_real_
    )
  })
  run <- grammar$run_study(s, 1)
  all <- grammar$assess_study(run)
  selected <- grammar$assess_study(run, function(measurements, attempts) {
    measurements$method_id == "first"
  })
  expect_equal(all$comparisons$n_used, c(2, 2))
  expect_equal(selected$comparisons$n_used, c(1, 1))
  expect_equal(selected$comparisons$n_failed, c(1, 1))
  expect_equal(nrow(selected$excluded), 2)
  empty <- grammar$assess_study(run, function(measurements, attempts) {
    rep(FALSE, nrow(measurements))
  })
  expect_equal(empty$comparisons$n_used, c(0, 0))
  expect_true(all(is.na(empty$comparisons$best)))
})

test_that("sequential and worker execution preserve results and caller RNG", {
  s <- make_grammar_study()
  withr::local_seed(402)
  before <- .Random.seed
  kind <- RNGkind()
  sequential <- grammar$run_study(s, 2)
  parallel <- grammar$run_study(s, 2, workers = 2)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)
  expect_equal(sequential$measurements, parallel$measurements)
})

test_that("measurements have separate random streams", {
  s <- make_grammar_study()
  s <- grammar$with_measure(s, "random", function(artifacts, context) {
    data.frame(value = runif(1))
  })
  first <- grammar$run_study(s, 2)
  s$measures <- s$measures[c("random", "estimate")]
  second <- grammar$run_study(s, 2)
  expect_equal(
    canonical_grammar_rows(first$measurements),
    canonical_grammar_rows(second$measurements)
  )
})

test_that("a failed comparison resumes from committed fits and transient inputs", {
  s <- make_grammar_study()
  marker <- tempfile()
  s$conditions$marker <- marker
  s$methods$first$fit <- function(data, context) {
    cat("fit\n", file = context$condition$marker, append = TRUE)
    list(data_mean = mean(data), draws = array(rnorm(40), c(20, 2, 1)))
  }
  s <- grammar$with_comparison(
    s,
    "interrupted",
    function(results, artifacts, context) {
      marker <- context$attempts$condition_marker[[1]]
      ready <- paste0(marker, ".ready")
      if (!file.exists(ready)) {
        file.create(ready)
        stop("comparison interrupted")
      }
      data.frame(value = mean(artifacts$first$draws))
    },
    needs = "draws"
  )
  s <- grammar$with_retention(s)
  path <- file.path(withr::local_tempdir(), "interrupted")
  expect_error(grammar$run_study(s, 1, path = path), "comparison interrupted")
  expect_length(readLines(marker), 1)
  resumed <- grammar$run_study(s, 1, path = path)
  expect_length(readLines(marker), 2) # only the second condition still needed a fit
  expect_equal(nrow(resumed$comparisons), 2)
  unlink(c(marker, paste0(marker, ".ready")))
})

test_that("planning rejects grouping fields the runtime cannot emit", {
  s <- make_grammar_study()
  s$generate <- function(condition, context) stop("must not generate")
  s$conditions$payload <- I(list(1:2, 3:4))
  s$conditions$matrix <- I(matrix(1:4, nrow = 2))
  for (id in names(s$methods)) {
    s$methods[[id]]$settings$vector <- 1:2
    s$methods[[id]]$settings$nested <- list(1)
    s$methods[[id]]$settings$scalar <- "ok"
  }
  for (by in c(
    "condition_payload",
    "condition_matrix",
    "method_vector",
    "method_nested"
  )) {
    invalid <- grammar$with_comparison(
      s,
      "invalid",
      function(results, artifacts, context) {
        data.frame()
      },
      by = by
    )
    expect_error(
      grammar$run_study(invalid, 1),
      "Unknown comparison grouping columns"
    )
  }
  valid <- grammar$with_comparison(
    s,
    "valid",
    function(results, artifacts, context) {
      data.frame()
    },
    by = c("condition_n", "method_scalar")
  )
  expect_equal(grammar$plan_study(valid, 1)$fits, 4)
})

test_that("retained reference artifacts are verified as bytes and corruption is rejected", {
  s <- make_grammar_study()
  s$methods <- s$methods["first"]
  s$methods$first$extract$reference <- function(fit, data, context) {
    shared <- new.env(parent = emptyenv())
    shared$draws <- fit$draws
    list(first = shared, second = shared)
  }
  s <- grammar$with_retention(s, c("draws", "truth", "reference"))
  path <- file.path(withr::local_tempdir(), "reference")
  first <- grammar$run_study(s, 1, path = path)
  s <- grammar$with_measure(
    s,
    "reference",
    function(artifacts, context) {
      stopifnot(identical(
        artifacts$reference$first,
        artifacts$reference$second
      ))
      data.frame(value = mean(artifacts$reference$first$draws))
    },
    needs = "reference"
  )
  reused <- grammar$run_study(s, 1, path = path)
  expect_identical(first$attempts, reused$attempts)
  expect_equal(nrow(reused$measurements), 4)
  file <- list.files(
    file.path(path, "fits"),
    recursive = TRUE,
    full.names = TRUE
  )[[1]]
  record <- readRDS(file)
  expect_type(record$payload, "raw")
  record$payload[[1]] <- as.raw(bitwXor(as.integer(record$payload[[1]]), 1L))
  saveRDS(record, file)
  expect_error(grammar$run_study(s, 1, path = path), "Checksum mismatch")
})

test_that("callback errors restore the caller's RNG state", {
  s <- make_grammar_study()
  s$generate <- function(condition, context) {
    runif(1)
    stop("generation interrupted")
  }
  withr::local_seed(492L)
  before <- .Random.seed
  kind <- RNGkind()
  expect_error(grammar$run_study(s, 1), "generation interrupted")
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)
})

test_that("the likelihood case selection implements the assessment policy contract", {
  example <- new.env(parent = grammar)
  sys.source(
    system.file("experimental", "likelihood-case.R", package = "bayesim"),
    example
  )
  s <- make_grammar_study()
  s <- grammar$with_measure(s, "performance", function(artifacts, context) {
    data.frame(
      value = 0.1,
      rmse_s = 0.5,
      comparable = TRUE,
      elpd_loo = -10,
      rhat = 1,
      ess_bulk = 1000,
      ess_tail = 1000,
      divergents = 0
    )
  })
  run <- grammar$run_study(s, 1)
  assessed <- grammar$assess_study(run, policy = example$likelihood_selection)
  expect_equal(nrow(assessed$measurements), 4)
  expect_true(all(assessed$measurements$measurement == "performance"))
  expect_equal(nrow(assessed$excluded), 4)
})
