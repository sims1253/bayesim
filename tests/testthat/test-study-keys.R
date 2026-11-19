# Keys, identities and what callbacks see. Invariant 2 of R/study-store.R:
# everything a callback can observe is in the key of what it produces.

one_method <- function(s, fit = function(data, context) list(v = data)) {
  s <- with_methods(s, m = study_method(fit))
  with_measure(s, "v", function(artifacts, context) {
    data.frame(value = 0, seen = artifacts$fit$v)
  }, needs = "fit")
}

# S15 ------------------------------------------------------------------------------

test_that("a summary of values that are all missing has no rows and the usual columns", {
  s <- study("na", data.frame(condition_id = "a"), function(condition, context) 1)
  s <- with_methods(s, m = study_method(function(data, context) list()))
  s <- with_measure(s, "v", function(artifacts, context) {
    data.frame(target = "t", value = NA_real_)
  }, needs = "fit")
  run <- run_study(s, 3)
  out <- summarise_study(run)
  expect_equal(nrow(out), 0)
  expect_named(out, c("condition_id", "method_id", "measurement", "target", "n", "estimate", "mcse"))
  aimed <- summarise_study(run, target_mcse = 0.1)
  expect_named(aimed, c(names(out), "n_needed", "done"))
})

# S17 ------------------------------------------------------------------------------

test_that("the serialised worker holds the run specification and nothing of its caller", {
  spec <- list(cap = list(), root = tempdir())
  bare <- length(serialize(bayesim:::.make_worker(spec), NULL))
  make <- function(spec) {
    big <- runif(131072) # 1 MB that the worker does not use
    bayesim:::.make_worker(spec)
  }
  expect_lt(length(serialize(make(spec), NULL)), bare + 10000)
})

# S7 -------------------------------------------------------------------------------

test_that("callbacks do not see the study name, so studies on one store share results", {
  seen <- withr::local_tempfile()
  make <- function(name) {
    s <- study(name, data.frame(condition_id = "a"), function(condition, context) {
      cat("generate", is.null(context$study), "\n", file = seen, append = TRUE)
      paste(context$study, "data")
    }, prepare = function(condition, context) {
      cat("prepare", is.null(context$study), "\n", file = seen, append = TRUE)
      1
    })
    s <- one_method(s)
    with_comparison(s, "c", function(results, artifacts, context) {
      cat("compare", is.null(context$study), "\n", file = seen, append = TRUE)
      data.frame(value = 0)
    }, needs = "fit")
  }
  path <- withr::local_tempfile()
  run_study(make("one"), 1, path = path)
  shared <- run_study(make("two"), 1, path = path)
  expect_equal(nrow(shared$executed), 0)
  expect_identical(shared$measurements$seen, run_study(make("two"), 1)$measurements$seen)
  expect_true(all(grepl("TRUE", readLines(seen))))
})

test_that("conditions that callbacks can tell apart have different keys", {
  make <- function(x) {
    s <- study("typed", data.frame(condition_id = "a", x = I(list(x))),
               function(condition, context) typeof(condition$x[[1]]))
    one_method(s)
  }
  make_column <- function(x) {
    s <- study("typed", data.frame(condition_id = "a", x = x),
               function(condition, context) paste(typeof(condition$x), class(condition$x)))
    one_method(s)
  }
  for (pair in list(list(1L, 1), list(TRUE, 1), list("a", factor("a")))) {
    path <- withr::local_tempfile()
    run_study(make_column(pair[[1]]), 1, path = path)
    second <- run_study(make_column(pair[[2]]), 1, path = path)
    expect_equal(nrow(second$executed), 1)
    expect_identical(second$measurements$seen, run_study(make_column(pair[[2]]), 1)$measurements$seen)
  }
  path <- withr::local_tempfile()
  run_study(make(1L), 1, path = path)
  expect_equal(nrow(run_study(make(1), 1, path = path)$executed), 1)
})

# S9 -------------------------------------------------------------------------------

test_that("the order of a measurement's needs is part of its key", {
  make <- function(needs) {
    s <- study("order", data.frame(condition_id = "a"), function(condition, context) 1)
    s <- with_methods(s, m = study_method(function(data, context) list()))
    s <- with_retention(s, "fit")
    with_measure(s, "first", function(artifacts, context) {
      data.frame(value = 0, first = names(artifacts)[1])
    }, needs = needs)
  }
  path <- withr::local_tempfile()
  expect_equal(run_study(make(c("fit", "data")), 1, path = path)$measurements$first, "fit")
  expect_equal(run_study(make(c("data", "fit")), 1, path = path)$measurements$first, "data")
})

test_that("comparisons with one callback and different names have their own rows and draws", {
  draw <- function(results, artifacts, context) data.frame(value = runif(1))
  s <- study("names", data.frame(condition_id = "a"), function(condition, context) 1)
  s <- with_methods(s, m = study_method(function(data, context) list()))
  s <- with_comparison(s, "one", draw, needs = "fit")
  s <- with_comparison(s, "two", draw, needs = "fit")
  run <- run_study(s, 1)
  expect_setequal(run$comparisons$comparison, c("one", "two"))
  expect_false(identical(
    run$comparisons$value[run$comparisons$comparison == "one"],
    run$comparisons$value[run$comparisons$comparison == "two"]
  ))
})

# S8 -------------------------------------------------------------------------------
#
# For each shape a factory builds a function from a value that the function
# reaches and a value that it does not. Changing the first must change the
# identity; changing the second must not.

identity_of <- function(f) bayesim:::.capture_function(f, "f")$id

expect_identity_tracks <- function(factory) {
  base <- identity_of(factory(1, 1))
  expect_identical(identity_of(factory(1, 1)), base)
  expect_false(identical(identity_of(factory(2, 1)), base))
  expect_identical(identity_of(factory(1, 2)), base)
}

test_that("a closure in a list counts with what it reaches", {
  expect_identity_tracks(function(a, b) {
    g <- local({
      a <- a
      b <- b
      function() a
    })
    lst <- list(g = g)
    function() lst$g()
  })
})

test_that("a closure in an environment counts with what it reaches", {
  expect_identity_tracks(function(a, b) {
    e <- new.env(parent = emptyenv())
    e$h <- local({
      a <- a
      b <- b
      function() a
    })
    function() e$h()
  })
})

test_that("a closure in an attribute counts with what it reaches", {
  expect_identity_tracks(function(a, b) {
    x <- structure(1, fn = local({
      a <- a
      b <- b
      function() a
    }))
    function() attr(x, "fn")()
  })
})

test_that("an object with methods in its own environment counts with its fields", {
  expect_identity_tracks(function(a, b) {
    unused <- b
    self <- new.env(parent = emptyenv())
    self$a <- a
    self$get <- function() self$a
    class(self) <- "counter"
    function() self$get()
  })
})

test_that("a list keeps its attributes in the identity", {
  expect_identity_tracks(function(a, b) {
    unused <- b
    x <- structure(list(1), meta = a)
    function() attr(x, "meta")
  })
})

test_that("a formula counts with the variables its environment gives it", {
  expect_identity_tracks(function(a, b) {
    fo <- local({
      k <- a
      z <- b
      y ~ x + k
    })
    function(data) stats::model.frame(fo, data)
  })
})

test_that("a package function counts by package and version, not by body", {
  versions <- c(stats = "1.0")
  local_mocked_bindings(.package_version = function(pkg) {
    if (pkg %in% names(versions)) versions[[pkg]] else "0"
  })
  direct <- identity_of(stats::median)
  free <- identity_of(function(x) median(x))
  qualified <- identity_of(function(x) stats::median(x))
  expect_identical(identity_of(stats::median), direct)
  versions[["stats"]] <- "2.0"
  expect_false(identical(identity_of(stats::median), direct))
  expect_false(identical(identity_of(function(x) median(x)), free))
  expect_false(identical(identity_of(function(x) stats::median(x)), qualified))
  # Two functions of one package with equal bodies are still two functions.
  expect_false(identical(identity_of(stats::median), identity_of(stats::quantile)))
})

test_that("primitives count by name", {
  make <- function(p) {
    ops <- list(p = p)
    function(x) ops$p(x)
  }
  expect_identical(identity_of(make(sum)), identity_of(make(sum)))
  expect_false(identical(identity_of(make(sum)), identity_of(make(prod))))
})

test_that("cycles end the walk and still count what they reach", {
  expect_identity_tracks(function(a, b) {
    unused <- b
    e <- new.env(parent = emptyenv())
    e$self <- e
    e$a <- a
    f <- function() e$self$a
    e$f <- f
    f
  })
})
