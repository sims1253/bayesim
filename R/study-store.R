# Content-addressed study store ---------------------------------------------
#
# Layout under the store directory:
#
#   bayesim-store                       marker, written once
#   fits/<xx>/<fit>-<measures>.rds      a fit and its measurements
#   fits/<xx>/<fit>.failed.rds          a fit that failed
#   artifacts/<xx>/<fit>.rds            artifacts kept after the fit
#   datasets/<xx>/<dataset>.rds         generated data, when retained
#   prepared/<xx>/<prepare>.rds         the value of prepare() for a condition
#   comparisons/<xx>/<key>.rds          one artifact-consuming comparison
#   claims/<key>/                       a job that some process runs now
#
# <xx> is the first two characters of the key. Each record holds a serialised
# payload and its checksum. Every file is written to a temporary name in its
# final directory and then renamed, so a reader sees the whole file or none of
# it. There is no index or manifest: the finished work is the set of file
# names, and a key changes when anything that determines the result changes.

.store_kinds <- c("fits", "artifacts", "datasets", "prepared", "comparisons")
.store_marker <- "bayesim study store, format 1"

.store_open <- function(path) {
  if (
    !is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)
  ) {
    cli::cli_abort("{.arg path} must be {.code NULL} or a directory path.")
  }
  if (file.exists(path) && !dir.exists(path)) {
    cli::cli_abort("{.arg path} {.path {path}} is a file, not a directory.")
  }
  marker <- file.path(path, "bayesim-store")
  if (file.exists(marker)) {
    if (!identical(readLines(marker, warn = FALSE), .store_marker)) {
      cli::cli_abort(
        "{.path {path}} uses a store format that this version of bayesim cannot read."
      )
    }
  } else {
    if (
      dir.exists(path) &&
        length(list.files(path, all.files = TRUE, no.. = TRUE))
    ) {
      cli::cli_abort(c(
        "{.path {path}} is not empty and is not a bayesim study store.",
        i = "Use a new directory."
      ))
    }
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    .write_atomic(
      .store_marker,
      marker,
      writer = function(x, f) writeLines(x, f)
    )
  }
  for (d in c(.store_kinds, "claims")) {
    dir.create(file.path(path, d), showWarnings = FALSE)
  }
  .store_handle(normalizePath(path, winslash = "/", mustWork = TRUE))
}

# The handle is a list of paths. Workers rebuild it from the root without IO.
.store_handle <- function(root) {
  dirs <- as.list(file.path(root, c(.store_kinds, "claims")))
  names(dirs) <- c(.store_kinds, "claims")
  c(list(root = root), dirs)
}

.store_file <- function(st, kind, key, suffix = ".rds") {
  file.path(st[[kind]], substr(key, 1L, 2L), paste0(key, suffix))
}

# File names present in the store, by kind, without the `.rds` suffix. One
# directory walk per kind.
.store_index <- function(st, kinds = .store_kinds) {
  out <- lapply(kinds, function(kind) {
    files <- list.files(st[[kind]], pattern = "\\.rds$", recursive = TRUE)
    sub("\\.rds$", "", basename(files))
  })
  names(out) <- kinds
  out
}

.write_atomic <- function(x, path, writer) {
  dir <- dirname(path)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  tmp <- file.path(
    dir,
    paste0(".tmp-", Sys.getpid(), "-", .next_tmp_id(), "-", basename(path))
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writer(x, tmp)
  if (!file.rename(tmp, path)) {
    cli::cli_abort("Could not write {.path {path}}.")
  }
  invisible(path)
}

# A per-process counter keeps temporary names unique without using the
# caller's random number stream.
.tmp_counter <- new.env(parent = emptyenv())
.next_tmp_id <- function() {
  .tmp_counter$n <- (.tmp_counter$n %||% 0) + 1
  format(.tmp_counter$n, scientific = FALSE)
}

.checksum <- function(payload) {
  digest::digest(payload, algo = "xxhash64", serialize = FALSE)
}

.write_payload <- function(value, path, compress = FALSE) {
  payload <- serialize(value, NULL, xdr = FALSE)
  .write_atomic(
    list(payload = payload, checksum = .checksum(payload)),
    path,
    writer = function(x, f) saveRDS(x, f, compress = compress)
  )
}

.read_payload <- function(path) {
  record <- tryCatch(readRDS(path), error = function(e) {
    cli::cli_abort("Unreadable study record {.path {path}}.", parent = e)
  })
  if (
    !is.list(record) ||
      !is.raw(record$payload) ||
      !identical(record$checksum, .checksum(record$payload))
  ) {
    cli::cli_abort(c(
      "Checksum mismatch: {.path {path}}.",
      i = "Delete the file to compute it again."
    ))
  }
  unserialize(record$payload)
}

# Read a record if it exists; NULL otherwise.
.read_if <- function(path) {
  if (file.exists(path)) .read_payload(path) else NULL
}

# Claims ---------------------------------------------------------------------
#
# A claim is a directory, because creating a directory is atomic and fails if
# it exists. A file in it holds the owner's host, process ID and start time. A
# claim is stale when its owner is a dead process on this host, or when it is
# older than `getOption("bayesim.claim_timeout")` seconds (default one day).

.host <- function() Sys.info()[["nodename"]]

.claim <- function(st, key) {
  dir <- file.path(st$claims, key)
  if (.claim_once(dir)) {
    return(TRUE)
  }
  if (.claim_stale(dir)) {
    unlink(dir, recursive = TRUE)
    return(.claim_once(dir))
  }
  FALSE
}

.claim_once <- function(dir) {
  if (!dir.create(dir, showWarnings = FALSE)) {
    return(FALSE)
  }
  writeLines(
    c(.host(), Sys.getpid(), format(unclass(Sys.time()), digits = 15)),
    file.path(dir, "owner")
  )
  TRUE
}

.claim_stale <- function(dir) {
  owner <- tryCatch(
    suppressWarnings(readLines(file.path(dir, "owner"), warn = FALSE)),
    error = function(e) character()
  )
  if (length(owner) < 3L) {
    # The owner may still be writing its details. Trust it for a minute.
    age <- as.numeric(Sys.time()) - as.numeric(file.mtime(dir))
    return(is.finite(age) && age > 60)
  }
  started <- as.numeric(owner[3])
  timeout <- getOption("bayesim.claim_timeout", 24 * 3600)
  if (is.finite(started) && as.numeric(Sys.time()) - started > timeout) {
    return(TRUE)
  }
  if (identical(owner[1], .host()) && .Platform$OS.type == "unix") {
    pid <- as.integer(owner[2])
    if (identical(pid, Sys.getpid())) {
      # Left by an earlier, interrupted call in this process.
      return(TRUE)
    }
    return(!isTRUE(tools::pskill(pid, signal = 0L)))
  }
  FALSE
}

.release <- function(st, key) {
  unlink(file.path(st$claims, key), recursive = TRUE)
}
