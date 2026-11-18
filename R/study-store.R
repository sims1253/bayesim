# Content-addressed study store ---------------------------------------------
#
# Invariants. The store, the jobs and the runs keep these. A change that
# breaks one is a bug, also when no test catches it.
#
# 1. The store is correct without claims. Claims save duplicate work; they
#    guard nothing. Two processes that run the same job at once must leave a
#    consistent store: every record complete, and every reference pointing to
#    the output of one and the same attempt.
# 2. A key names one result. Everything a callback can observe is in the
#    key, and nothing else is.
# 3. A run reports what it used. It decides once, per job, which record it
#    uses, collects exactly those records, and does not list the store again
#    at the end.
# 4. Retention trades space for time, never for correctness. An input that
#    was discarded and is needed again is computed again.
# 5. A job sees the same inputs on every backend and in every order. No job
#    can observe what another job did to a shared object.
#
# How the store keeps invariant 1:
#
# * Every file is written once. A writer writes a temporary file and links
#   it to its final name. The link fails when the name exists, so of two
#   writers exactly one wins, and the other uses the file that won.
# * Every fit attempt has an ID. The artifacts of an attempt are files named
#   by that ID, and the record of the attempt holds it. A comparison key
#   names the attempts that the comparison read. A record therefore never
#   points to the artifacts of another attempt, and the loser of a race
#   deletes only its own files.
# * A record is published after everything it promises: its retained
#   artifacts and, when the data are retained, its dataset.
#
# Layout under the store directory:
#
#   bayesim-store                          marker, written once
#   fits/<xx>/<fit>-<measures>.rds         a successful attempt and its
#                                          measurements
#   fits/<xx>/<fit>.failed.<attempt>.rds   a failed attempt; never deleted
#   artifacts/<xx>/<fit>.<attempt>.<h>.rds one retained artifact of an attempt
#   transient/<xx>/<fit>.<attempt>.<h>.rds one artifact that only
#                                          comparisons need; the comparison
#                                          job deletes it
#   datasets/<xx>/<dataset>.rds            generated data, when retained
#   prepared/<xx>/<prepare>.rds            the value of prepare() for a
#                                          condition
#   comparisons/<xx>/<key>.rds             one comparison result
#   claims/<key>/<token>                   a job that some process runs now
#
# <xx> is the first two characters of the key, and <h> the first eight of
# the artifact name's digest. Each record holds a serialised payload and its
# checksum. There is no index or manifest: the finished work is the set of
# file names, and a key changes when anything that determines the result
# changes.

.store_kinds <- c(
  "fits", "artifacts", "transient", "datasets", "prepared", "comparisons"
)

# Store trace ----------------------------------------------------------------
#
# Every read, write, deletion, directory walk, probe, claim and release goes
# through a function that reports itself to .trace_op(). The tests count
# these operations to check how much work a run does without timing it. When
# no trace is on, the cost is one test of a logical.

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

.store_marker <- "bayesim study store, format 2"

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

# Fit records -------------------------------------------------------------------
#
# A fit key has 40 characters. The record of a successful attempt is
# `<fit>-<measures>`, and that of a failed attempt `<fit>.failed.<attempt>`.

.fit_key_of <- function(names) substr(names, 1L, 40L)
.is_success_record <- function(names) substr(names, 41L, 41L) == "-"
.is_failure_record <- function(names) substr(names, 41L, 48L) == ".failed."

.failure_record <- function(fit_key, attempt_id) {
  paste0(fit_key, ".failed.", attempt_id)
}

# The failure records of one fit, from one listing of its shard.
.fit_failures <- function(st, fit_key) {
  dir <- file.path(st$fits, substr(fit_key, 1L, 2L))
  .trace_op("probe", dir)
  files <- list.files(
    dir,
    pattern = paste0("^", fit_key, "\\.failed\\..*\\.rds$")
  )
  sub("\\.rds$", "", files)
}

# The file of one artifact of one attempt, in `kind` "artifacts" or
# "transient".
.artifact_file <- function(st, kind, fit_key, attempt_id, name) {
  .store_file(
    st,
    kind,
    paste(fit_key, attempt_id, substr(.hash_text("artifact", name), 1L, 8L), sep = ".")
  )
}

# Attempt IDs sort by the time they were made: milliseconds since the epoch,
# then the process nonce and a counter. The smallest ID of a fit is its
# canonical attempt.
.attempt_id <- function() {
  paste0(
    sprintf("%013.0f", floor(as.numeric(Sys.time()) * 1000)),
    .process_nonce(),
    sprintf("%09d", as.integer(.next_tmp_id()))
  )
}

# Writing ----------------------------------------------------------------------

# Write `x` to `path` once. The writer writes a temporary file in the final
# directory, then links it to `path`. A link fails when `path` exists, so of
# several processes that write one path, exactly one wins. Returns TRUE when
# this call put the file in place and FALSE when the file was already there.
#
# On a file system without hard links, the temporary file is renamed
# instead. A rename replaces an existing file, so there the last writer wins.
.write_atomic <- function(x, path, writer) {
  dir <- dirname(path)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  tmp <- file.path(dir, .tmp_name(basename(path)))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  .trace_op("write", path)
  writer(x, tmp)
  if (suppressWarnings(file.link(tmp, path))) {
    return(invisible(TRUE))
  }
  if (file.exists(path)) {
    return(invisible(FALSE))
  }
  if (!file.rename(tmp, path)) {
    cli::cli_abort("Could not write {.path {path}}.")
  }
  invisible(TRUE)
}

# Delete store files that this process wrote and no longer uses.
.store_delete <- function(paths) {
  for (p in paths) .trace_op("delete", p)
  unlink(paths)
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

# A record is an RDS file of a list that holds the serialised payload and
# its checksum. Without compression its bytes are those of serialize(), so
# both sides skip the connection that saveRDS() and readRDS() open.
# Returns TRUE when this call wrote the record and FALSE when it existed.
.write_raw <- function(payload, path, compress = FALSE) {
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

.write_payload <- function(value, path, compress = FALSE) {
  .write_raw(serialize(value, NULL, xdr = FALSE), path, compress)
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

# The checked payload bytes of a record. `size` lets a caller that reads
# many records stat them in one call. A size that does not match the file
# makes the record fail its check; the file is then read again with a fresh
# size. Records are written once, so this happens only when a user replaced
# a file.
.read_raw <- function(path, size = NULL) {
  record <- tryCatch(.read_record(path, size), error = function(e) e)
  if (!is.null(size) && !.record_ok(record)) {
    return(.read_raw(path))
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
  record$payload
}

.read_payload <- function(path, size = NULL) {
  unserialize(.read_raw(path, size))
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

# Publish `value` at `path` unless a record is there. Returns the payload
# bytes of the record that is in the store afterwards: this call's, or the
# one that another process wrote first.
.publish_raw <- function(payload, path, compress = FALSE) {
  if (.write_raw(payload, path, compress)) payload else .read_raw(path)
}

# Claims ---------------------------------------------------------------------
#
# A claim says that some process runs a job now, so that other processes
# can do other work. By invariant 1, a claim guards nothing: a job that two
# processes run at once costs only duplicate work.
#
# A claim is a directory `claims/<key>` that holds one empty file, the
# token. The token names its owner: `<host>_<pid>_<ms>_<nonce>`. To claim,
# a process makes a directory with its token under a temporary name and
# renames it to the claim's name. The rename fails when the claim exists
# and holds a token. To take over a stale claim, a process renames the old
# token to its own. A rename of one file succeeds once, so of several
# processes that take over one claim, exactly one wins.
#
# A claim is stale when this process made it and no longer holds it, when
# its owner is a dead process on this host, or when it is older than
# `getOption("bayesim.claim_timeout")` seconds (default one day). A process
# releases only a claim that holds its own token.

.claims_held <- new.env(parent = emptyenv())

.host <- function() {
  .tmp_counter$host <- .tmp_counter$host %||%
    gsub("[^A-Za-z0-9.-]", "-", Sys.info()[["nodename"]])
  .tmp_counter$host
}

.claim_token <- function() {
  paste(
    .host(),
    Sys.getpid(),
    sprintf("%.0f", floor(as.numeric(Sys.time()) * 1000)),
    paste0(.process_nonce(), "-", .next_tmp_id()),
    sep = "_"
  )
}

# Claim the job `key`. Returns the token, or NULL when another process holds
# the claim.
.claim <- function(st, key) {
  dir <- file.path(st$claims, key)
  token <- .claim_token()
  stage <- file.path(st$claims, paste0(".tmp-", token))
  dir.create(stage, showWarnings = FALSE)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  file.create(file.path(stage, token))
  won <- function() {
    .claims_held[[token]] <- key
    .trace_op("claim", dir)
    token
  }
  for (try in 1:5) {
    # A rename onto a missing or empty directory succeeds.
    if (suppressWarnings(file.rename(stage, dir))) {
      return(won())
    }
    owners <- list.files(dir)
    if (!length(owners)) {
      # Released meanwhile.
      next
    }
    old <- owners[[1L]]
    if (!.claim_stale(old)) {
      return(NULL)
    }
    if (suppressWarnings(file.rename(file.path(dir, old), file.path(dir, token)))) {
      return(won())
    }
  }
  NULL
}

.claim_stale <- function(token) {
  owner <- strsplit(token, "_", fixed = TRUE)[[1L]]
  if (length(owner) != 4L) {
    # Not a token of this format. Its owner cannot be known.
    return(TRUE)
  }
  pid <- suppressWarnings(as.integer(owner[[2L]]))
  here <- identical(owner[[1L]], .host())
  if (here && identical(pid, Sys.getpid())) {
    return(is.null(.claims_held[[token]]))
  }
  started <- as.numeric(owner[[3L]]) / 1000
  timeout <- getOption("bayesim.claim_timeout", 24 * 3600)
  if (!is.finite(started) || as.numeric(Sys.time()) - started > timeout) {
    return(TRUE)
  }
  if (here && .Platform$OS.type == "unix" && !is.na(pid)) {
    return(!isTRUE(tools::pskill(pid, signal = 0L)))
  }
  FALSE
}

# Whether a live claim on `key` exists.
.claim_busy <- function(st, key) {
  owners <- list.files(file.path(st$claims, key))
  length(owners) > 0L && !.claim_stale(owners[[1L]])
}

.release <- function(st, key, token) {
  if (is.null(token)) {
    return(invisible(FALSE))
  }
  if (!is.null(.claims_held[[token]])) {
    rm(list = token, envir = .claims_held)
  }
  dir <- file.path(st$claims, key)
  .trace_op("release", dir)
  mine <- suppressWarnings(file.remove(file.path(dir, token)))
  if (mine) {
    # Fails when another process has claimed the job meanwhile.
    suppressWarnings(file.remove(dir))
  }
  invisible(mine)
}
