# test-metric-registry.R
# Tests for metric registry functions: register_metric, get_metric, clear_registry, etc.

library(bayesim)

describe("Metric Registry", {
  describe("register_metric()", {
    it("registers a metric by name", {
      clear_registry()
      m <- rmse_metric(name = "test_rmse")
      expect_silent(register_metric(m))
      expect_true("test_rmse" %in% list_metrics())
    })

    it("rejects non-S7 objects", {
      clear_registry()
      expect_error(
        register_metric(list(name = "not_s7")),
        "metric must be an S7 object"
      )
    })

    it("rejects non-Metric S7 objects", {
      clear_registry()
      f <- MockFitter()
      expect_error(
        register_metric(f),
        "metric must inherit from the Metric S7 class"
      )
    })

    it("rejects metrics with empty name", {
      clear_registry()
      # Create a metric with empty name via S7
      EmptyMetric <- S7::new_class(
        "EmptyMetric",
        parent = Metric,
        properties = list(
          name = S7::new_property(S7::class_character, default = "")
        )
      )
      m <- EmptyMetric()
      # S7 may return NA or empty for default - either way it should error
      expect_error(register_metric(m))
    })

    it("rejects duplicate registration by default", {
      clear_registry()
      m1 <- rmse_metric(name = "dup_test")
      m2 <- bias_metric(name = "dup_test")
      register_metric(m1)
      expect_error(
        register_metric(m2),
        "already registered"
      )
    })

    it("allows overwrite when explicitly requested", {
      clear_registry()
      m1 <- rmse_metric(name = "overwrite_test")
      m2 <- bias_metric(name = "overwrite_test")
      register_metric(m1)
      expect_silent(register_metric(m2, overwrite = TRUE))
      retrieved <- get_metric("overwrite_test")
      # S7 class names include package prefix
      expect_true(inherits(retrieved, "BiasMetric") || grepl("BiasMetric", class(retrieved)[1]))
    })
  })

  describe("get_metric()", {
    it("retrieves a registered metric by name", {
      clear_registry()
      m <- rmse_metric(name = "get_test")
      register_metric(m)
      retrieved <- get_metric("get_test")
      expect_true(inherits(retrieved, "RmseMetric") || grepl("RmseMetric", class(retrieved)[1]))
      expect_equal(retrieved@name, "get_test")
    })

    it("returns NULL for unregistered metric", {
      clear_registry()
      expect_null(get_metric("nonexistent_metric"))
    })

    it("rejects non-character name", {
      clear_registry()
      expect_error(get_metric(123), "name must be a single character string")
    })

    it("rejects vector name", {
      clear_registry()
      expect_error(get_metric(c("a", "b")), "name must be a single character string")
    })
  })

  describe("list_metrics()", {
    it("returns empty vector when registry is empty", {
      clear_registry()
      expect_equal(list_metrics(), character())
    })

    it("returns names of all registered metrics", {
      clear_registry()
      register_metric(rmse_metric(name = "lm1"))
      register_metric(bias_metric(name = "lm2"))
      register_metric(coverage_metric(name = "lm3"))
      metrics <- list_metrics()
      expect_true(setequal(metrics, c("lm1", "lm2", "lm3")))
    })
  })

  describe("unregister_metric()", {
    it("removes a metric from the registry", {
      clear_registry()
      register_metric(rmse_metric(name = "unreg_test"))
      expect_true("unreg_test" %in% list_metrics())
      unregister_metric("unreg_test")
      expect_false("unreg_test" %in% list_metrics())
    })

    it("warns when unregistering nonexistent metric", {
      clear_registry()
      expect_warning(
        unregister_metric("nonexistent"),
        "not registered"
      )
    })

    it("rejects non-character name", {
      clear_registry()
      expect_error(unregister_metric(123), "name must be a single character string")
    })
  })

  describe("clear_registry()", {
    it("removes all metrics", {
      clear_registry()
      register_metric(rmse_metric(name = "clr1"))
      register_metric(bias_metric(name = "clr2"))
      expect_true(length(list_metrics()) > 0)
      clear_registry()
      expect_equal(list_metrics(), character())
    })

    it("is safe to call on empty registry", {
      clear_registry()
      expect_silent(clear_registry())
    })
  })

  describe("resolve_metrics_from_registry()", {
    it("resolves character names to Metric objects", {
      clear_registry()
      register_metric(rmse_metric(name = "rr1"))
      register_metric(bias_metric(name = "rr2"))
      resolved <- resolve_metrics_from_registry(c("rr1", "rr2"))
      expect_equal(length(resolved), 2)
      expect_true(inherits(resolved[[1]], "RmseMetric") || grepl("RmseMetric", class(resolved[[1]])[1]))
      expect_true(inherits(resolved[[2]], "BiasMetric") || grepl("BiasMetric", class(resolved[[2]])[1]))
    })

    it("returns empty list for NULL input", {
      clear_registry()
      expect_equal(resolve_metrics_from_registry(NULL), list())
    })

    it("passes through list of Metric objects", {
      clear_registry()
      metrics <- list(rmse_metric(name = "passthrough1"), bias_metric(name = "passthrough2"))
      resolved <- resolve_metrics_from_registry(metrics)
      expect_equal(length(resolved), 2)
    })

    it("errors on unregistered character name", {
      clear_registry()
      expect_error(
        resolve_metrics_from_registry("nonexistent"),
        "not found in registry"
      )
    })

    it("rejects non-Metric objects in list", {
      clear_registry()
      expect_error(
        resolve_metrics_from_registry(list("not_an_object")),
        "not an S7 Metric object"
      )
    })

    it("rejects non-character non-list input", {
      clear_registry()
      expect_error(
        resolve_metrics_from_registry(123),
        "must be a character vector or list"
      )
    })
  })
})
