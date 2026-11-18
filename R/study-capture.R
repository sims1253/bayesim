# Closure capture and identity -----------------------------------------------
#
# A study's user functions (generator, methods, measurements, comparisons)
# must run on
# worker processes that do not share the controller's global environment. At
# plan time each function is rebuilt as a self-contained closure: its free
# variables are looked up the way R would look them up, and the values found
# outside package namespaces are copied into a fresh environment whose parent
# is the base environment. Functions found in package namespaces are bound by
# reference; they serialise as a namespace reference and resolve on the worker
# when the package is installed there.
#
# The identity of a function is a digest of .stable(): its formals and body
# (without source references or comments), and everything it can reach
# through its free variables, nested closures, environments and formulas
# included. Package functions count by package, version and name.

.capture_limit <- function() {
  getOption("bayesim.capture_limit", 10 * 1024^2)
}

# Returns list(fn = self-contained function, id = identity digest).
.capture_function <- function(f, label) {
  state <- new.env(parent = emptyenv())
  state$seen <- list()
  fn <- .capture_rec(f, label, label, state)$fn
  list(fn = fn, id = .value_id(f))
}

.capture_rec <- function(f, label, path, state) {
  if (!is.function(f)) {
    cli::cli_abort("{label} must be a function.", call = NULL)
  }
  if (is.primitive(f)) {
    return(list(fn = f))
  }
  env <- environment(f)
  if (.is_shared_env(env)) {
    # Package function passed directly: it travels as a namespace reference.
    return(list(fn = f))
  }
  for (s in state$seen) {
    if (identical(s$original, f)) {
      return(list(fn = s$fn))
    }
  }
  crate_env <- new.env(parent = baseenv())
  crated <- utils::removeSource(f)
  environment(crated) <- crate_env
  slot <- length(state$seen) + 1L
  state$seen[[slot]] <- list(original = f, fn = crated)

  names <- sort(unique(unlist(
    codetools::findGlobals(f, merge = FALSE),
    use.names = FALSE
  )))
  for (name in names) {
    home <- .find_binding(name, env)
    if (is.null(home) || .is_base_env(home)) {
      # Unresolved names are often column names used through NSE; base
      # bindings exist on every worker.
      next
    }
    value <- get(name, envir = home, inherits = FALSE)
    if (.is_shared_env(home)) {
      assign(name, value, envir = crate_env)
      next
    }
    where <- paste0(path, "$", name)
    if (is.function(value)) {
      inner <- .capture_rec(value, label, where, state)
      assign(name, inner$fn, envir = crate_env)
    } else {
      .check_capturable(value, label, name)
      assign(name, value, envir = crate_env)
    }
  }
  list(fn = crated)
}

.fn_shape <- function(f) {
  f <- utils::removeSource(f)
  list(formals = as.list(formals(f)), body = body(f))
}

.find_binding <- function(name, env) {
  while (!identical(env, emptyenv())) {
    if (exists(name, envir = env, inherits = FALSE)) {
      return(env)
    }
    env <- parent.env(env)
  }
  NULL
}

.is_base_env <- function(env) {
  identical(env, baseenv()) || identical(env, .BaseNamespaceEnv)
}

# Namespaces, their import environments and attached package environments are
# available on workers by name and are never copied.
.is_shared_env <- function(env) {
  isNamespace(env) ||
    .is_base_env(env) ||
    identical(env, emptyenv()) ||
    startsWith(environmentName(env), "imports:") ||
    startsWith(environmentName(env), "package:")
}

.env_label <- function(env) {
  if (isNamespace(env)) getNamespaceName(env) else environmentName(env)
}

# Identity of a captured value: a digest of everything a function can reach
# through it. See .stable().
.value_id <- function(value) {
  .digest(list("value", .stable(value)))
}

# A form of `x` that digests the same in every session and changes when
# anything reachable from `x` changes:
#
# * Atomic values and lists carry their attributes, classed lists included.
# * A closure counts by its formals and body, the packages its body names
#   with `::`, and the values of its free variables, looked up the way R
#   looks them up. Values in other bindings of its environment do not count.
# * A package function counts by package, version and name. Its body does
#   not count, because the version names it.
# * A primitive counts by name.
# * A formula counts by its text and the values that its environment gives
#   to the names in it.
# * An environment counts by all its bindings and its parents. Namespaces,
#   package environments and the global environment count by name only.
# * Closures and environments that the walk has already entered count by
#   their place in the walk, so cycles end.
.stable <- function(x, state = .stable_state()) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.function(x)) {
    return(.stable_function(x, state))
  }
  if (is.environment(x)) {
    return(.stable_env(x, state))
  }
  if (inherits(x, "formula")) {
    return(.stable_formula(x, state))
  }
  attrs <- .stable_attrs(attributes(x), state)
  if (isS4(x)) {
    return(list("S4", typeof(x), attrs))
  }
  if (is.list(x)) {
    y <- unclass(x)
    items <- lapply(seq_along(y), function(i) .stable(y[[i]], state))
    return(list("list", typeof(x), items, attrs))
  }
  if (is.null(attrs)) {
    return(x)
  }
  attributes(x) <- NULL
  list("attributes", x, attrs)
}

.stable_state <- function() {
  state <- new.env(parent = emptyenv())
  state$seen <- list()
  state
}

.stable_attrs <- function(attrs, state) {
  attrs$srcref <- NULL
  attrs$srcfile <- NULL
  attrs$wholeSrcref <- NULL
  if (!length(attrs)) {
    return(NULL)
  }
  attrs <- attrs[order(names(attrs))]
  lapply(attrs, .stable, state = state)
}

# The place of `x` in the walk if the walk has entered it, or NULL after
# entering it now.
.stable_enter <- function(x, state) {
  for (i in seq_along(state$seen)) {
    if (identical(state$seen[[i]], x)) {
      return(list("seen", i))
    }
  }
  state$seen[[length(state$seen) + 1L]] <- x
  NULL
}

.stable_function <- function(f, state) {
  if (is.primitive(f)) {
    return(list("primitive", .primitive_name(f)))
  }
  env <- environment(f)
  if (.is_shared_env(env)) {
    return(.stable_package_function(f, env))
  }
  seen <- .stable_enter(f, state)
  if (!is.null(seen)) {
    return(seen)
  }
  names <- sort(unique(unlist(
    codetools::findGlobals(f, merge = FALSE),
    use.names = FALSE
  )))
  free <- list()
  for (name in names) {
    home <- .find_binding(name, env)
    if (is.null(home) || .is_base_env(home)) {
      next
    }
    free[[name]] <- if (.is_shared_env(home)) {
      # An import environment holds other packages' functions; such a
      # function counts with the package that defines it.
      value <- .get_binding(name, home)
      owner <- if (is.function(value) && !is.primitive(value)) environment(value)
      pkg <- .env_package(if (!is.null(owner) && .is_shared_env(owner)) owner else home)
      list("package", pkg, .package_version(pkg), name)
    } else {
      .stable(.get_binding(name, home), state)
    }
  }
  pkgs <- sort(unique(.qualified_packages(body(f))))
  list(
    "closure",
    .fn_shape(f),
    free,
    lapply(stats::setNames(nm = pkgs), .package_version),
    .stable_attrs(attributes(f), state)
  )
}

.stable_package_function <- function(f, env) {
  pkg <- .env_package(env)
  name <- .namespace_name(f, env)
  list(
    "package",
    pkg,
    .package_version(pkg),
    name %||% .fn_shape(f)
  )
}

.stable_env <- function(e, state) {
  if (identical(e, globalenv())) {
    return(list("environment", "global"))
  }
  if (.is_shared_env(e)) {
    return(list("environment", .env_label(e)))
  }
  seen <- .stable_enter(e, state)
  if (!is.null(seen)) {
    return(seen)
  }
  names <- sort(ls(e, all.names = TRUE))
  values <- lapply(names, function(name) {
    if (bindingIsActive(name, e)) "active binding" else .stable(.get_binding(name, e), state)
  })
  list(
    "environment",
    names,
    values,
    .stable_env(parent.env(e), state),
    .stable_attrs(attributes(e), state)
  )
}

.stable_formula <- function(x, state) {
  env <- environment(x)
  values <- list()
  if (is.environment(env)) {
    for (name in sort(unique(all.names(x)))) {
      home <- .find_binding(name, env)
      if (is.null(home) || .is_shared_env(home)) {
        next
      }
      values[[name]] <- .stable(.get_binding(name, home), state)
    }
  }
  attrs <- attributes(x)
  attrs$.Environment <- NULL
  text <- x
  attributes(text) <- NULL
  list("formula", deparse(text), values, .stable_attrs(attrs, state))
}

# A binding's value. A binding can hold the missing argument, which get()
# cannot return.
.get_binding <- function(name, env) {
  tryCatch(
    get(name, envir = env, inherits = FALSE),
    error = function(e) list("unreadable binding", conditionMessage(e))
  )
}

.primitive_name <- function(f) {
  sub('^\\.Primitive\\("(.*)"\\)$', "\\1", deparse(f)[[1L]])
}

.env_package <- function(env) {
  if (.is_base_env(env)) {
    return("base")
  }
  sub("^(package|imports):", "", .env_label(env))
}

# The name under which a namespace holds `f`, or NULL.
.namespace_name <- function(f, env) {
  if (!isNamespace(env)) {
    return(NULL)
  }
  for (name in sort(ls(env, all.names = TRUE))) {
    if (!bindingIsActive(name, env) &&
        identical(get(name, envir = env, inherits = FALSE), f)) {
      return(name)
    }
  }
  NULL
}

# Package names that an expression uses with `::` or `:::`.
.qualified_packages <- function(expr) {
  if (!is.call(expr) && !is.pairlist(expr)) {
    return(character())
  }
  out <- character()
  if (is.call(expr) && length(expr) == 3L &&
      (identical(expr[[1L]], quote(`::`)) || identical(expr[[1L]], quote(`:::`)))) {
    out <- as.character(expr[[2L]])
  }
  # Index the parts: a part can be the empty symbol of a missing argument,
  # which cannot be bound to a loop variable and used.
  parts <- as.list(expr)
  for (i in seq_along(parts)) {
    if (is.call(parts[[i]]) || (is.pairlist(parts[[i]]) && length(parts[[i]]))) {
      out <- c(out, .qualified_packages(parts[[i]]))
    }
  }
  out
}

# The installed version of a package, or "" when it is not installed. The
# result is kept for the life of the process.
.package_versions <- new.env(parent = emptyenv())
.package_version <- function(pkg) {
  if (is.null(.package_versions[[pkg]])) {
    .package_versions[[pkg]] <- tryCatch(
      as.character(utils::packageVersion(pkg)),
      error = function(e) ""
    )
  }
  .package_versions[[pkg]]
}

.check_capturable <- function(value, label, name) {
  bad <- .find_unserialisable(value, new.env(parent = emptyenv()))
  if (!is.null(bad)) {
    cli::cli_abort(
      c(
        "{label} captures {.var {name}}, which cannot be sent to a worker.",
        i = "{.var {name}} contains {bad}, which does not survive serialisation.",
        i = "Create it inside the function instead of capturing it."
      ),
      call = NULL
    )
  }
  size <- length(serialize(value, NULL))
  if (size > .capture_limit()) {
    cli::cli_abort(
      c(
        "{label} captures {.var {name}}, which is {prettyNum(size, big.mark = ',')} bytes when serialised.",
        i = "Each job sends captured values to its worker. The limit is {prettyNum(.capture_limit(), big.mark = ',')} bytes.",
        i = "Read large inputs inside the function, or raise {.code options(bayesim.capture_limit = )}."
      ),
      call = NULL
    )
  }
  invisible(NULL)
}

# Depth-first search for external pointers. Returns a description or NULL.
.find_unserialisable <- function(x, visited) {
  type <- typeof(x)
  if (type == "externalptr") {
    return("an external pointer")
  }
  if (inherits(x, "connection")) {
    return("a connection")
  }
  if (type == "environment") {
    if (.is_shared_env(x) || identical(x, globalenv())) {
      return(NULL)
    }
    tag <- format(x)
    if (!is.null(visited[[tag]])) {
      return(NULL)
    }
    visited[[tag]] <- TRUE
    for (nm in ls(x, all.names = TRUE)) {
      if (bindingIsActive(nm, x)) {
        next
      }
      hit <- .find_unserialisable(get(nm, envir = x, inherits = FALSE), visited)
      if (!is.null(hit)) {
        return(hit)
      }
    }
    return(.find_unserialisable(parent.env(x), visited))
  }
  if (type == "closure") {
    return(.find_unserialisable(environment(x), visited))
  }
  if (is.list(x)) {
    for (el in x) {
      hit <- .find_unserialisable(el, visited)
      if (!is.null(hit)) {
        return(hit)
      }
    }
  }
  for (a in attributes(x)) {
    hit <- .find_unserialisable(a, visited)
    if (!is.null(hit)) {
      return(hit)
    }
  }
  NULL
}

# Capture every function of a study once and compute the identities that
# keys are built from. The result is what workers receive.
.sg_capture <- function(study) {
  label <- function(...) paste0(...)
  gen <- .capture_function(study$generate, "The generator")
  prep <- if (!is.null(study$prepare)) {
    .capture_function(study$prepare, "The prepare function")
  }
  methods <- lapply(names(study$methods), function(id) {
    m <- study$methods[[id]]
    fit <- .capture_function(m$fit, label("Method \"", id, "\" (fit)"))
    extract <- lapply(names(m$extract), function(nm) {
      .capture_function(
        m$extract[[nm]],
        label("Method \"", id, "\" (extractor \"", nm, "\")")
      )
    })
    names(extract) <- names(m$extract)
    .check_capturable(m$settings, label("Method \"", id, "\""), "settings")
    list(
      fit = fit$fn,
      extract = lapply(extract, `[[`, "fn"),
      settings = m$settings,
      fit_id = .digest(list(fit$id, .value_id(m$settings), m$version)),
      extract_id = .digest(list(
        lapply(extract, `[[`, "id")[sort(names(extract))],
        m$version
      ))
    )
  })
  names(methods) <- names(study$methods)
  measures <- lapply(names(study$measures), function(nm) {
    x <- study$measures[[nm]]
    f <- .capture_function(x$compute, label("Measurement \"", nm, "\""))
    list(
      compute = f$fn,
      needs = x$needs,
      summary = x$summary %||% "mean",
      # A callback receives its artifacts in the order of `needs`, so the
      # order is part of the identity.
      id = .digest(list(f$id, x$needs, x$version))
    )
  })
  names(measures) <- names(study$measures)
  comparisons <- lapply(names(study$comparisons), function(nm) {
    x <- study$comparisons[[nm]]
    f <- .capture_function(x$compute, label("Comparison \"", nm, "\""))
    list(
      compute = f$fn,
      by = x$by,
      needs = x$needs,
      id = .digest(list(f$id, x$by, x$needs, x$version))
    )
  })
  names(comparisons) <- names(study$comparisons)
  list(
    name = study$name,
    generate = gen$fn,
    prepare = prep$fn,
    generate_id = .digest(list(gen$id, study$version)),
    prepare_id = .digest(list(prep$id, study$version)),
    methods = methods,
    measures = measures,
    comparisons = comparisons,
    retention = study$retention,
    compress = study$compress %||% "gzip",
    required = .sg_required(study)
  )
}
