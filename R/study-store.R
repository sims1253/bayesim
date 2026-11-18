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

# Store trace ----------------------------------------------------------------
#
# Every read, write, directory walk and claim goes through one of four
# functions, and each reports itself to .trace_op(). The tests count these
# operations to check how much work a run does without timing it. When no
# trace is on, the cost is one test of a logical.

.store_trace <- new.env(parent = emptyenv())
.store_trace$on <- FALSE

.trace_op <- function(op, path) {
  if (.store_trace$on) {
    n <- .store_trace$n + 1L
    .store_trace$op[[n]] <- op
    .store_trace$path[[n]] <- path
    .store_trace$n <- n
  }
}

# Evaluate `expr` and return its value with the store operations that this
# process made meanwhile, as a data frame of `op` and `path`. Jobs that a
# backend runs in another process are not counted.
.with_store_trace <- function(expr) {
  .store_trace$n <- 0L
  .store_trace$op <- character()
  .store_trace$path <- character()
  .store_trace$on <- TRUE
  on.exit(.store_trace$on <- FALSE)
  value <- expr
  list(
    value = value,
    ops = data.frame(
      op = .store_trace$op,
      path = .store_trace$path,
      stringsAsFactors = FALSE
    )
  )
}
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
    # An earlier open that stopped before the marker was in place leaves
    # only entries that the store itself makes. They are not foreign.
    entries <- list.files(path, all.files = TRUE, no.. = TRUE)
    own <- entries %in% c(.store_kinds, "claims", "bayesim-store") |
      startsWith(entries, ".tmp-")
    if (dir.exists(path) && !all(own)) {
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
    .trace_op("walk", st[[kind]])
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
  tmp <- file.path(dir, .tmp_name(basename(path)))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  .trace_op("write", path)
  writer(x, tmp)
  if (!file.rename(tmp, path)) {
    cli::cli_abort("Could not write {.path {path}}.")
  }
  invisible(path)
}

# A temporary name is unique across hosts and processes: it holds the host,
# the process ID, a nonce for this process and a counter. Stores on shared
# file systems see many hosts, and process IDs repeat across hosts and over
# time.
.tmp_counter <- new.env(parent = emptyenv())
.tmp_name <- function(base) {
  paste0(
    ".tmp-", .host(), "-", Sys.getpid(), "-", .process_nonce(), "-",
    .next_tmp_id(), "-", base
  )
}

# The nonce is a digest of the host, the process ID, the time to the
# microsecond and the name of the session's temporary directory. It draws no
# random numbers, so it does not move the caller's stream. A forked child
# has another process ID and makes its own nonce.
.process_nonce <- function() {
  pid <- Sys.getpid()
  if (!identical(.tmp_counter$pid, pid)) {
    .tmp_counter$nonce <- substr(
      .vdigest("sha1")(paste(
        .host(),
        pid,
        format(Sys.time(), "%Y-%m-%d %H:%M:%OS6"),
        proc.time()[["elapsed"]],
        basename(tempdir()),
        sep = "\r"
      )),
      1L,
      16L
    )
    .tmp_counter$n <- 0L
    .tmp_counter$pid <- pid
  }
  .tmp_counter$nonce
}

.next_tmp_id <- function() {
  .process_nonce()
  .tmp_counter$n <- .tmp_counter$n + 1L
  as.character(.tmp_counter$n)
}

.checksum <- function(payload) .vdigest("xxhash64")(payload)

# A record is an RDS file. Without compression its bytes are those of
# serialize(), so both sides skip the connection that saveRDS() and readRDS()
# open.
.write_payload <- function(value, path, compress = FALSE) {
  payload <- serialize(value, NULL, xdr = FALSE)
  record <- list(payload = payload, checksum = .checksum(payload))
  .write_atomic(
    record,
    path,
    writer = if (isFALSE(compress)) {
      function(x, f) writeBin(serialize(x, NULL), f)
    } else {
      function(x, f) saveRDS(x, f, compress = compress)
    }
  )
}

.read_record <- function(path, size = NULL) {
  .trace_op("read", path)
  size <- size %||% file.size(path)
  bytes <- if (is.finite(size) && size > 2) readBin(path, "raw", size)
  if (length(bytes) == size && identical(bytes[1:2], as.raw(c(0x58, 0x0a)))) {
    unserialize(bytes)
  } else {
    readRDS(path)
  }
}

# `size` lets a caller that reads many records stat them in one call. A
# size taken before the file changed reads too few or too many bytes; the
# record then fails its check, and the file is read again with a fresh size.
.read_payload <- function(path, size = NULL) {
  record <- tryCatch(.read_record(path, size), error = function(e) e)
  if (!is.null(size) && !.record_ok(record)) {
    return(.read_payload(path))
  }
  if (inherits(record, "error")) {
    cli::cli_abort("Unreadable study record {.path {path}}.", parent = record)
  }
  if (!.record_ok(record)) {
    cli::cli_abort(c(
      "Checksum mismatch: {.path {path}}.",
      i = "Delete the file to compute it again."
    ))
  }
  unserialize(record$payload)
}

.record_ok <- function(record) {
  is.list(record) &&
    is.raw(record$payload) &&
    identical(record$checksum, .checksum(record$payload))
}

# Whether store files exist. Every probe of the store for a record goes
# through here, so that the trace counts it.
.store_has <- function(paths) {
  for (p in paths) .trace_op("probe", p)
  file.exists(paths)
}

# Read a record if it exists; NULL otherwise.
.read_if <- function(path) {
  if (.store_has(path)) .read_payload(path) else NULL
}

# Claims ---------------------------------------------------------------------
#
# A claim is a directory, because creating a directory is atomic and fails if
# it exists. The name of the one empty file in it holds the owner's host,
# process ID and start time. A claim is stale when its owner is a dead process
# on this host, or when it is older than
# `getOption("bayesim.claim_timeout")` seconds (default one day).

.host <- function() {
  .tmp_counter$host <- .tmp_counter$host %||%
    gsub("[^A-Za-z0-9.-]", "-", Sys.info()[["nodename"]])
  .tmp_counter$host
}

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
  .trace_op("claim", dir)
  file.create(file.path(
    dir,
    paste(.host(), Sys.getpid(), sprintf("%.3f", unclass(Sys.time())), sep = "_")
  ))
  TRUE
}

.claim_stale <- function(dir) {
  owner <- strsplit(list.files(dir)[1L], "_", fixed = TRUE)[[1L]]
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
