# Derived seeds and keys -----------------------------------------------------
#
# Nothing here draws from a random stream in sequence. Each seed and key is a
# digest of the values that define it. A result therefore does not depend on
# job order, the number of workers, the backend, or which other jobs exist.

.digest <- function(x) digest::digest(x, algo = "sha1")

# Vectorised digests, made once per process.
.digesters <- new.env(parent = emptyenv())
.vdigest <- function(algo) {
  fn <- .digesters[[algo]]
  if (is.null(fn)) {
    vd <- digest::getVDigest(algo)
    fn <- function(x) vd(x, serialize = FALSE)
    .digesters[[algo]] <- fn
  }
  fn
}

# One sha1 per element. The arguments are recycled and joined with a
# separator that cannot occur in a hexadecimal key or a decimal number.
.hash_text <- function(...) {
  .vdigest("sha1")(paste(..., sep = "\r"))
}

# The key of each condition row. All fields count, including condition_id.
# The key hashes the row exactly as the functions of the study receive it:
# with column order, types and attributes. 1L and 1 are different conditions,
# because a function can tell them apart.
.condition_keys <- function(rows) {
  vapply(rows, .value_id, character(1), USE.NAMES = FALSE)
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
# and state, also when `expr` fails. The first word of `.Random.seed` names
# the generator, so assigning the state also sets the kind. Only a caller
# without a state needs RNGkind() to restore its kind. Under
# `normal.kind = "Box-Muller"` the generator caches a second normal draw
# outside `.Random.seed`, so a restored caller loses that one draw.
.with_stream <- function(seed, expr) {
  env <- globalenv()
  old_seed <- env[[".Random.seed"]]
  old_kind <- if (is.null(old_seed)) RNGkind()
  on.exit({
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = env)
    } else {
      RNGkind(old_kind[1], old_kind[2], old_kind[3])
      if (exists(".Random.seed", envir = env, inherits = FALSE)) {
        rm(".Random.seed", envir = env)
      }
    }
  })
  assign(".Random.seed", c(10407L, seed), envir = env)
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
