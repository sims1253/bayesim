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
# The same walk yields the function's identity: a digest of its formals, its
# body (without source references or comments), and the identities of every
# captured value. Package functions are identified by package and name only,
# so changes inside installed packages do not change identities.

.capture_limit <- function() {
  getOption("bayesim.capture_limit", 10 * 1024^2)
}

# Returns list(fn = self-contained function, id = identity digest).
.capture_function <- function(f, label) {
  state <- new.env(parent = emptyenv())
  state$seen <- list()
  .capture_rec(f, label, label, state)
}

.capture_rec <- function(f, label, path, state) {
  if (!is.function(f)) {
    cli::cli_abort("{label} must be a function.", call = NULL)
  }
  if (is.primitive(f)) {
    return(list(fn = f, id = .digest(list("primitive", deparse(f)))))
  }
  env <- environment(f)
  if (.is_shared_env(env)) {
    # Package function passed directly: it travels as a namespace reference.
    return(list(fn = f, id = .digest(list("package", .fn_shape(f)))))
  }
  for (s in state$seen) {
    if (identical(s$original, f)) {
      return(list(fn = s$fn, id = s$id %||% .digest(list("cycle", path))))
    }
  }
  crate_env <- new.env(parent = baseenv())
  crated <- utils::removeSource(f)
  environment(crated) <- crate_env
  slot <- length(state$seen) + 1L
  state$seen[[slot]] <- list(original = f, fn = crated, id = NULL)

  names <- sort(unique(unlist(
    codetools::findGlobals(f, merge = FALSE),
    use.names = FALSE
  )))
  ids <- list()
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
      ids[[name]] <- list("package", .env_label(home), name)
      next
    }
    where <- paste0(path, "$", name)
    if (is.function(value)) {
      inner <- .capture_rec(value, label, where, state)
      assign(name, inner$fn, envir = crate_env)
      ids[[name]] <- inner$id
    } else {
      .check_capturable(value, label, name)
      assign(name, value, envir = crate_env)
      ids[[name]] <- .value_id(value)
    }
  }
  id <- .digest(list("closure", .fn_shape(crated), ids))
  state$seen[[slot]]$id <- id
  list(fn = crated, id = id)
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

# Identity of a captured value. Formulas count by their text, because their
# environment changes from session to session.
.value_id <- function(value) {
  .digest(list("value", .stable_value(value)))
}

.stable_value <- function(value) {
  if (inherits(value, "formula")) {
    return(list("formula", deparse(value)))
  }
  if (is.function(value)) {
    return(list("function", .fn_shape(value)))
  }
  if (is.list(value) && !is.data.frame(value) && is.null(attr(value, "class"))) {
    return(lapply(value, .stable_value))
  }
  value
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
