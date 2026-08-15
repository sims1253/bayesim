# mori shared-memory integration for the model bank.
#
# When bayesim runs on LOCAL daemons, the model bank is an identical read-only
# payload shipped to every worker via mirai::everywhere(). Without mori each
# daemon deserializes its own full copy; mori::share() writes the bank into a
# single OS shared-memory region that every daemon maps zero-copy instead.
#
# mori requires shared physical RAM, so it is correct ONLY for local daemons.
# A remote daemon would receive a serialized shared-memory *reference* (e.g.
# "/mori_<pid>_<id>") to a region on the controller's machine that it cannot
# map, silently corrupting get_model_bank() on that worker. Locality is
# therefore decided conservatively: anything that is not provably local
# (loopback transport) is treated as remote and the bank is passed through
# unchanged.

# Classify a mirai daemon URL as local (same-machine) or not.
#
# Machine-local transports — abstract://, ipc://, inproc:// — can only connect
# processes on the same host. Network transports (tcp://, tls+tcp://, ws://,
# wss://) are local only when bound to a loopback address. Wildcard binds
# (0.0.0.0) and any non-loopback host are treated as remote: this is
# conservative (the worst case is a missed optimization, never corruption).
#
# `url` is the value returned by mirai::status()$daemons.
#
# @param url Scalar character, e.g. "abstract://abc123" or "tcp://127.0.0.1:5555".
# @return TRUE if the transport is provably local, FALSE otherwise.
# @keywords internal
# @noRd
.is_local_daemon_url <- function(url) {
  if (is.na(url) || !nzchar(url)) {
    return(FALSE)
  }
  # Machine-local transports.
  if (
    startsWith(url, "abstract://") ||
      startsWith(url, "ipc://") ||
      startsWith(url, "inproc://")
  ) {
    return(TRUE)
  }
  # Network transports: strip scheme, then strip trailing :port (and IPv6
  # brackets) to recover the host, and accept only loopback addresses.
  rest <- sub("^[a-zA-Z+]+://", "", url)
  host <- sub(":[0-9]+$", "", rest)
  host <- sub("^\\[(.+)\\]$", "\\1", host)
  host %in% c("127.0.0.1", "localhost", "::1", "0:0:0:0:0:0:0:1")
}

# Report whether the currently-set mirai daemons are all on the local machine.
#
# Queries mirai::status() for the daemon URL and classifies it. Returns FALSE
# (not "share-eligible") when no daemons are set, when the URL cannot be read,
# or when the transport is not provably local — so callers can gate mori
# sharing on this single predicate without a separate sequential-mode check.
#
# @return TRUE if mirai daemons are set and use a local transport; FALSE
#   otherwise (including the no-daemons / sequential case).
# @keywords internal
# @noRd
.are_daemons_local <- function() {
  if (!isTRUE(mirai::daemons_set())) {
    return(FALSE)
  }
  url <- tryCatch(
    mirai::status()$daemons,
    error = function(e) NA_character_
  )
  .is_local_daemon_url(url)
}

# Share the model bank across local daemons using mori shared memory.
#
# Returns mori::share(bank) when local daemons are detected (so each daemon
# zero-copy maps the same physical pages instead of receiving a serialized
# copy), and returns `bank` unchanged otherwise. Remote-daemon and sequential
# runs are never affected.
#
# Lifetime: the caller must keep the returned object reachable for as long as
# daemons need to map it (i.e. for the duration of the run). mori frees the
# shared-memory region once no R object — in any process — references it. In
# execute_tasks() the return value is held in a local variable for the whole
# run, satisfying this contract.
#
# mori is optional because current releases require a newer R than bayesim's
# supported minimum. When it is absent the ordinary serialized mirai dispatch
# remains correct. If shared-region creation fails at runtime, the same safe
# fallback is used with a warning.
#
# @param bank A model bank (named list) or NULL.
# @return The (possibly shared) bank, or NULL if `bank` is NULL. The original
#   `bank` is never mutated.
# @keywords internal
# @noRd
.maybe_share_bank <- function(bank) {
  if (is.null(bank)) {
    return(NULL)
  }
  if (!.are_daemons_local()) {
    return(bank)
  }
  if (!requireNamespace("mori", quietly = TRUE)) {
    return(bank)
  }
  tryCatch(
    {
      shared <- mori::share(bank)
      cli::cli_alert_info(
        "Model bank shared in local memory via mori (daemons map it zero-copy)."
      )
      shared
    },
    error = function(e) {
      cli::cli_warn(c(
        "mori sharing unavailable; falling back to serialized dispatch.",
        i = conditionMessage(e)
      ))
      bank
    }
  )
}
