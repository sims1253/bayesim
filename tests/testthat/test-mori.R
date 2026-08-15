# test-mori.R — mori shared-memory model-bank integration

describe(".is_local_daemon_url", {
  it("classifies machine-local transports as local", {
    expect_true(bayesim:::.is_local_daemon_url("abstract://abc123"))
    expect_true(bayesim:::.is_local_daemon_url("ipc:///tmp/sock"))
    expect_true(bayesim:::.is_local_daemon_url("inproc://x"))
  })

  it("classifies loopback tcp/ws/tls as local", {
    expect_true(bayesim:::.is_local_daemon_url("tcp://127.0.0.1:5555"))
    expect_true(bayesim:::.is_local_daemon_url("tcp://localhost:5555"))
    expect_true(bayesim:::.is_local_daemon_url("tcp://[::1]:5555"))
    expect_true(bayesim:::.is_local_daemon_url("tls+tcp://127.0.0.1:5555"))
    expect_true(bayesim:::.is_local_daemon_url("ws://127.0.0.1:5555"))
    expect_true(bayesim:::.is_local_daemon_url("wss://[::1]:5555"))
  })

  it("classifies wildcard / non-loopback / hostname as NOT local (conservative)", {
    # 0.0.0.0 binds all interfaces — daemons could be remote, so treat as remote.
    expect_false(bayesim:::.is_local_daemon_url("tcp://0.0.0.0:5555"))
    expect_false(bayesim:::.is_local_daemon_url("tcp://10.0.0.5:5555"))
    expect_false(bayesim:::.is_local_daemon_url("tcp://192.168.1.10:5555"))
    expect_false(bayesim:::.is_local_daemon_url("tcp://myhost:5555"))
    expect_false(bayesim:::.is_local_daemon_url("ws://example.com:5555"))
    expect_false(bayesim:::.is_local_daemon_url("tls+tcp://10.0.0.5:5555"))
  })

  it("handles missing / empty / malformed input", {
    expect_false(bayesim:::.is_local_daemon_url(NA_character_))
    expect_false(bayesim:::.is_local_daemon_url(""))
  })
})

describe(".are_daemons_local", {
  it("returns FALSE when no daemons are set", {
    mirai::daemons(0)
    on.exit(mirai::daemons(0), add = TRUE)
    expect_false(bayesim:::.are_daemons_local())
  })

  it("returns TRUE for default local daemons (abstract socket)", {
    mirai::daemons(2)
    on.exit(mirai::daemons(0), add = TRUE)
    expect_true(bayesim:::.are_daemons_local())
  })

  it("returns FALSE for a non-loopback daemon URL (via mocked status)", {
    # Force the locality predicate to see a remote URL without binding a socket.
    local_mocked_bindings(
      status = function() {
        list(daemons = "tcp://10.0.0.5:5555", connections = 4L)
      },
      daemons_set = function() TRUE,
      .package = "mirai"
    )
    expect_false(bayesim:::.are_daemons_local())
  })
})

describe(".maybe_share_bank", {
  it("returns NULL for a NULL bank", {
    expect_null(bayesim:::.maybe_share_bank(NULL))
  })

  it("passes the bank through unchanged when no daemons are set", {
    mirai::daemons(0)
    on.exit(mirai::daemons(0), add = TRUE)
    bank <- list(spec1 = list(values = seq_len(100)))
    expect_identical(bayesim:::.maybe_share_bank(bank), bank)
  })

  it("passes the bank through unchanged for remote daemons", {
    local_mocked_bindings(
      status = function() {
        list(daemons = "tcp://10.0.0.5:5555", connections = 4L)
      },
      daemons_set = function() TRUE,
      .package = "mirai"
    )
    bank <- list(spec1 = list(values = seq_len(100)))
    # No shared region should be created for remote daemons.
    expect_identical(bayesim:::.maybe_share_bank(bank), bank)
  })

  it("passes the bank through unchanged when the optional package is unavailable", {
    local_mocked_bindings(
      .are_daemons_local = function() TRUE,
      .package = "bayesim"
    )
    local_mocked_bindings(
      requireNamespace = function(package, quietly = FALSE) FALSE,
      .package = "base"
    )
    bank <- list(spec1 = list(values = seq_len(100)))
    expect_identical(bayesim:::.maybe_share_bank(bank), bank)
  })

  it("shares the bank when local daemons are set", {
    skip_if_not_installed("mori")
    mirai::daemons(2)
    on.exit(mirai::daemons(0), add = TRUE)
    bank <- list(spec1 = list(values = seq_len(10000)))
    out <- bayesim:::.maybe_share_bank(bank)
    # A real shared-memory region is created (name is non-NULL) ...
    expect_false(is.null(mori::shared_name(out)))
    # ... and the content survives the round-trip into shared memory.
    expect_identical(out$spec1$values, bank$spec1$values)
  })

  it("end-to-end: a daemon reads the shared bank correctly", {
    skip_if_not_installed("mori")
    mirai::daemons(2)
    on.exit(mirai::daemons(0), add = TRUE)
    bank <- list(spec1 = list(values = seq_len(5000)))
    shared <- bayesim:::.maybe_share_bank(bank)
    expect_false(is.null(mori::shared_name(shared)))
    # Ship exactly as execute_tasks() does.
    mirai::everywhere(
      options(bayesim.model_bank = mb),
      .args = list(mb = shared)
    )
    # A daemon maps the shared region and reads element 5000. (No cleanup
    # everywhere() needed: mirai::daemons(0) in the on.exit destroys the
    # workers, and the shared region is GC'd when `shared` goes out of scope.)
    res <- mirai::mirai(getOption("bayesim.model_bank")$spec1$values[5000])[]
    expect_equal(res, 5000)
  })
})
