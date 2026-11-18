# Derived seeds and keys -----------------------------------------------------
#
# Nothing here draws from a random stream in sequence. Each seed and key is a
# digest of the values that define it. A result therefore does not depend on
# job order, the number of workers, the backend, or which other jobs exist.

.digest <- function(x) digest::digest(x, algo = "sha1")

.vdigest <- function(algo) {
  fn <- digest::getVDigest(algo)
  function(x) fn(x, serialize = FALSE)
}

# One sha1 per element. The arguments are recycled and joined with a
# separator that cannot occur in a hexadecimal key or a decimal number.
.hash_text <- function(...) {
  .vdigest("sha1")(paste(..., sep = "\r"))
}

# The key of each condition row. All fields count, including condition_id.
# Whole numbers stored as integers and as doubles give the same key, and
# factors count as their labels.
.condition_keys <- function(conditions) {
  vapply(
    seq_len(nrow(conditions)),
    function(i) {
      row <- lapply(as.list(conditions[i, , drop = FALSE]), function(x) {
        x <- unclass(x)
        if (is.factor(x)) {
          x <- as.character(x)
        }
        if (is.integer(x) || is.logical(x)) {
          storage.mode(x) <- "double"
        }
        x
      })
      .digest(row[order(names(row))])
    },
    character(1)
  )
}

# The L'Ecuyer-CMRG state for a stream name. Each of the six seed words takes
# 28 bits of a SHA-256 digest, so both halves stay below their moduli and are
# never all zero.
.stream_seed <- function(stream) {
  h <- .vdigest("sha256")(paste(stream, collapse = "\r"))
  starts <- (0:5) * 7L + 1L
  as.integer(strtoi(substring(h, starts, starts + 6L), 16L) + 1L)
}

# Evaluate `expr` on the stream `seed`, then restore the caller's generator
# and state, also when `expr` fails.
.with_stream <- function(seed, expr) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
  on.exit({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  assign(".Random.seed", c(10407L, seed), envir = globalenv())
  expr
}

# Call a user function on its own stream. `context$seed` receives the first
# word of the stream, for functions that pass a seed to other software.
.invoke <- function(stream, fn, ...) {
  seed <- .stream_seed(stream)
  args <- list(...)
  if ("context" %in% names(args)) {
    args$context$seed <- seed[[1L]]
  }
  .with_stream(seed, do.call(fn, args, quote = TRUE))
}
