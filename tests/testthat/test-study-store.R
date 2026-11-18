# The study store: publishing, claims, attempts and what a run reads. The
# invariants these tests check are written at the top of R/study-store.R.

local_store <- function(env = parent.frame()) {
  bayesim:::.store_open(withr::local_tempdir(.local_envir = env))
}

fit_name <- function(i = 1L) strrep(sprintf("%x", i %% 16L), 40L)

# S5 -------------------------------------------------------------------------------

test_that("temporary names hold the host, the process and a nonce, and draw no random numbers", {
  st <- local_store()
  seen <- character()
  writer <- function(x, f) {
    seen <<- c(seen, basename(f))
    writeLines(x, f)
  }
  set.seed(1)
  before <- .Random.seed
  bayesim:::.write_atomic("a", file.path(st$root, "one"), writer)
  bayesim:::.write_atomic("b", file.path(st$root, "two"), writer)
  expect_identical(.Random.seed, before)
  host <- gsub(".", "\\.", bayesim:::.host(), fixed = TRUE)
  pattern <- paste0("^\\.tmp-", host, "-", Sys.getpid(), "-([0-9a-f]{16})-[0-9]+-(one|two)$")
  expect_match(seen, pattern)
  nonce <- sub(pattern, "\\1", seen[[1]])
  # Another process gets another nonce.
  skip_if_not_installed("callr")
  other <- callr::r(function() {
    set.seed(1)
    before <- .Random.seed
    nonce <- bayesim:::.process_nonce()
    list(nonce = nonce, same_seed = identical(.Random.seed, before))
  })
  expect_true(other$same_seed)
  expect_false(identical(other$nonce, nonce))
})

# S3 -------------------------------------------------------------------------------

test_that("opening a store again after an interrupted first open succeeds", {
  root <- withr::local_tempdir()
  # What the first open leaves when it stops before the marker is in place.
  dir.create(file.path(root, "fits"))
  dir.create(file.path(root, "claims"))
  file.create(file.path(root, ".tmp-host-1-0123456789abcdef-1-bayesim-store"))
  st <- bayesim:::.store_open(root)
  expect_true(file.exists(file.path(st$root, "bayesim-store")))
  foreign <- withr::local_tempdir()
  file.create(file.path(foreign, "notes.txt"))
  expect_error(bayesim:::.store_open(foreign), "not empty")
})

# S6 -------------------------------------------------------------------------------

test_that("a size taken before a record changed does not read as corruption", {
  st <- local_store()
  path <- bayesim:::.store_file(st, "datasets", fit_name())
  bayesim:::.write_payload(list(x = 1), path)
  size <- file.size(path)
  unlink(path)
  bayesim:::.write_payload(list(x = 1:1000), path)
  expect_identical(bayesim:::.read_payload(path, size)$x, 1:1000)
})
