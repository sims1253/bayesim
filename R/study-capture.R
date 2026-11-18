# Capture and identity --------------------------------------------------------
#
# A study's callbacks run on worker processes that do not share the
# controller's global environment. At planning, .capture() turns each
# callback, each method's settings and each condition row into a transport
# form. The identity that keys are built from is a digest of that form, so
# what is hashed is what travels: a binding counts in the identity exactly
# when a worker receives it.
#
# The transport form of a value:
#
# * A closure whose environment is not a package's gets a copy of its
#   environment chain. Each copy holds only the bindings that the captured
#   code reads there, and its parent is the copy of the original parent.
#   Closures that share an environment share its copy. The copy of the
#   global environment has the base environment as its parent.
# * A name that the code reads from a package (a namespace, its imports or
#   an attached package) is bound by reference in the last copy of the
#   chain. It travels as a namespace reference and counts by package,
#   version and name. A name that base defines is not copied.
# * An environment reached as a value is copied with all its bindings.
# * A formula's environment is copied with the bindings that the names in
#   the formula read there.
# * Lists, expressions, calls and attributes are rewritten part by part.
#   Attributes keep their order.
#
# .capture() refuses, and names the callback and the binding, anything that
# it cannot carry exactly: an active binding, a promise of a name or a call
# that is not forced,
# a missing argument, the `...` of an enclosing function, the global
# environment as a value, an external pointer, a weak reference, a
# connection, and a form larger than `getOption("bayesim.capture_limit")`.
#
# The names that a closure reads are its free variables, as codetools finds
# them, and every name that its code can read before it assigns it. A
# replacement such as `x$a <- 1` reads `x` from the enclosing environment
# before it makes a local copy, and `x <- x + 1` reads `x` before it
# assigns it. Copying a binding that the code never reads is safe; missing
# one is not.

.capture_limit <- function() {
  getOption("bayesim.capture_limit", 10 * 1024^2)
}

# The transport form of a callback and its identity.
.capture_function <- function(f, label) {
  if (!is.function(f)) {
    cli::cli_abort("{label} must be a function.", call = NULL)
  }
  out <- .capture(f, label)
  list(fn = out$value, id = out$id)
}

# The transport form of any value and its identity. `name` names the value
# in messages when the value is not a function.
.capture <- function(x, label, name = "") {
  state <- new.env(parent = emptyenv())
  state$label <- label
  state$envs <- list()
  state$closures <- list()
  state$filled <- new.env(parent = emptyenv())
  state$bindings <- list()
  value <- .capture_value(x, state, name)
  .capture_check_size(value, state, name)
  list(value = value, id = .value_id(value))
}

# `problem` and `hint` are cli templates; `{name}` in them is the binding.
.capture_refuse <- function(state, path, problem, hint, name = "", code = "") {
  label <- state$label
  what <- if (nzchar(path)) path else "a value"
  cli::cli_abort(
    c(
      "{label} captures {.var {what}}, which cannot be sent to a worker.",
      i = problem,
      i = hint
    ),
    call = NULL
  )
}

.capture_path <- function(path, name) {
  if (nzchar(path)) paste0(path, "$", name) else name
}

.capture_value <- function(x, state, path) {
  if (is.null(x)) {
    return(x)
  }
  type <- typeof(x)
  if (type == "externalptr") {
    .capture_refuse(state, path, "It holds an external pointer, which does not survive serialisation.",
                    "Create it inside the function instead of capturing it.")
  }
  if (type == "weakref") {
    .capture_refuse(state, path, "It holds a weak reference, which does not survive serialisation.",
                    "Create it inside the function instead of capturing it.")
  }
  if (inherits(x, "connection")) {
    .capture_refuse(state, path, "It is a connection, which does not survive serialisation.",
                    "Open the connection inside the function instead of capturing it.")
  }
  if (type == "closure") {
    return(.capture_closure(x, state, path))
  }
  if (type %in% c("builtin", "special")) {
    return(x)
  }
  if (type == "environment") {
    return(.capture_env_value(x, state, path))
  }
  if (type %in% c("list", "expression", "language")) {
    x <- .capture_parts(x, state, path)
  }
  .capture_attributes(x, state, path)
}

# Replace the parts of a list, expression or call that capture changes.
.capture_parts <- function(x, state, path) {
  parts <- as.list(x)
  keys <- names(parts)
  changed <- FALSE
  for (i in seq_along(parts)) {
    if (is.null(parts[[i]]) || is.symbol(parts[[i]])) {
      next
    }
    key <- if (!is.null(keys) && nzchar(keys[[i]])) keys[[i]] else paste0("[[", i, "]]")
    old <- parts[[i]]
    new <- .capture_value(old, state, .capture_path(path, key))
    if (!identical(new, old)) {
      parts[[i]] <- new
      changed <- TRUE
    }
  }
  if (!changed) {
    return(x)
  }
  if (is.call(x)) {
    out <- as.call(parts)
  } else {
    out <- parts
    if (is.expression(x)) out <- as.expression(out)
  }
  attributes(out) <- attributes(x)
  out
}

# Replace the attributes that capture changes, in place, so that their
# order stays.
.capture_attributes <- function(x, state, path) {
  attrs <- attributes(x)
  for (nm in names(attrs)) {
    if (nm %in% c("srcref", "srcfile", "wholeSrcref")) {
      next
    }
    old <- attrs[[nm]]
    where <- paste0("attr(", if (nzchar(path)) path else "value", ", \"", nm, "\")")
    new <- if (nm == ".Environment" && is.environment(old) && is.language(x)) {
      .capture_formula_env(x, old, state, where)
    } else {
      .capture_value(old, state, where)
    }
    if (!identical(new, old)) {
      attr(x, nm) <- new
    }
  }
  x
}

.capture_closure <- function(f, state, path) {
  env <- environment(f)
  if (.is_shared_env(env)) {
    # A package function travels as a namespace reference.
    return(f)
  }
  for (s in state$closures) {
    if (identical(s$original, f)) {
      return(s$copy)
    }
  }
  copy <- utils::removeSource(f)
  environment(copy) <- .capture_crate(env, state)$copy
  slot <- length(state$closures) + 1L
  state$closures[[slot]] <- list(original = f, copy = copy)
  if (.uses_free_dots(f) && !is.null(.find_binding("...", env))) {
    .capture_refuse(
      state, .capture_path(path, "..."),
      "It is the `...` of an enclosing function, whose parts are promises.",
      "Collect them with {.code list(...)} in the enclosing function and capture the list."
    )
  }
  for (name in .closure_names(f)) {
    .capture_name(name, env, state, path)
  }
  copy <- .capture_attributes(copy, state, path)
  state$closures[[slot]]$copy <- copy
  copy
}

# The copy of `env`. Copies are made once per environment, and registered
# before they are filled, so that cycles end.
.capture_crate <- function(env, state) {
  for (s in state$envs) {
    if (identical(s$original, env)) {
      return(s)
    }
  }
  copy <- new.env(parent = emptyenv())
  slot <- length(state$envs) + 1L
  entry <- list(original = env, copy = copy, slot = slot, whole = FALSE)
  state$envs[[slot]] <- entry
  if (identical(env, globalenv())) {
    parent.env(copy) <- baseenv()
    return(entry)
  }
  parent <- parent.env(env)
  parent.env(copy) <- if (.is_shared_env(parent)) parent else .capture_crate(parent, state)$copy
  attrs <- attributes(env)
  for (nm in names(attrs)) {
    attr(copy, nm) <- .capture_value(attrs[[nm]], state, paste0("attr(environment, \"", nm, "\")"))
  }
  if (isS4(env)) {
    copy <- asS4(copy)
    entry$copy <- copy
    state$envs[[slot]] <- entry
  }
  entry
}

# An environment reached as a value travels with all its bindings.
.capture_env_value <- function(e, state, path) {
  if (.is_shared_env(e)) {
    return(e)
  }
  if (identical(e, globalenv())) {
    .capture_refuse(state, path, "It is the global environment, which a worker does not share.",
                    "Capture the values that the function needs instead.")
  }
  entry <- .capture_crate(e, state)
  if (!entry$whole) {
    entry$whole <- TRUE
    state$envs[[entry$slot]] <- entry
    names <- .radix_sort(ls(e, all.names = TRUE, sorted = FALSE))
    for (name in names) {
      .capture_fill(entry, name, e, state, .capture_path(path, name))
    }
  }
  entry$copy
}

# A formula's environment, with the bindings that the formula's names read.
.capture_formula_env <- function(x, env, state, path) {
  if (.is_shared_env(env)) {
    return(env)
  }
  for (name in .radix_sort(unique(all.names(x)))) {
    .capture_name(name, env, state, path)
  }
  .capture_crate(env, state)$copy
}

# Copy the binding that `name` reads from `env`, as R looks it up, into
# the copy of the environment that holds it. A binding past the end of the
# copied chain (a package, or the search path above the global environment)
# goes to the last copy of the chain.
.capture_name <- function(name, env, state, path) {
  inside <- TRUE
  last <- NULL
  home <- env
  repeat {
    if (identical(home, emptyenv())) {
      # Not bound anywhere: a column name used through non-standard
      # evaluation, or an error that the function raises when it runs.
      return(invisible())
    }
    if (inside && .is_shared_env(home)) {
      inside <- FALSE
    }
    if (inside) {
      last <- home
    }
    if (exists(name, envir = home, inherits = FALSE)) {
      break
    }
    if (identical(home, globalenv())) {
      inside <- FALSE
    }
    home <- parent.env(home)
  }
  if (.is_base_env(home)) {
    return(invisible())
  }
  entry <- .capture_crate(if (inside) home else last, state)
  .capture_fill(entry, name, home, state, .capture_path(path, name))
}

.capture_fill <- function(entry, name, home, state, path) {
  key <- paste(entry$slot, name, sep = "\r")
  if (!is.null(state$filled[[key]])) {
    return(invisible())
  }
  state$filled[[key]] <- TRUE
  value <- if (.is_shared_env(home)) {
    # Attached packages bind lazily loaded data; reading forces it.
    get(name, envir = home, inherits = FALSE)
  } else {
    .capture_binding(name, home, state, path)
  }
  value <- .capture_value(value, state, path)
  assign(name, value, envir = entry$copy)
  if (!.is_shared_env(home)) {
    state$bindings[[path]] <- value
  }
  invisible()
}

# The value of a binding that capture copies, or a refusal.
.capture_binding <- function(name, env, state, path) {
  if (name == "..." || grepl("^[.][.][0-9]+$", name)) {
    .capture_refuse(state, path, "It is the `...` of a function, whose parts are promises.",
                    "Collect them with {.code list(...)} and capture the list.")
  }
  if (bindingIsActive(name, env)) {
    .capture_refuse(state, path, "It is an active binding, which computes its value each time it is read.",
                    "Assign its value to an ordinary variable and capture that.")
  }
  # A promise whose code is a value, as from `make(1)` or do.call(), gives
  # that value whatever happens before it is forced. A promise of a name or
  # a call can give another value later, or fail, so capture refuses it.
  code <- if (rlang::env_binding_are_lazy(env, name)) {
    do.call(substitute, list(as.name(name), env))
  }
  # do.call() and callr pass values as `quote(<value>)`, which reads nothing.
  quoted <- is.call(code) &&
    (identical(code[[1L]], quote(quote)) || identical(code[[1L]], quote(base::quote)))
  if (is.language(code) && !quoted) {
    .capture_refuse(state, path, "It is a promise that is not forced yet, so it reads {.code {code}} when the callback first runs.",
                    "Call {.code force({name})} in the function that makes the callback.",
                    name = name, code = paste(deparse(code), collapse = " "))
  }
  missing <- new.env(parent = emptyenv())
  value <- tryCatch(get(name, envir = env, inherits = FALSE), error = function(e) missing)
  if (identical(value, missing)) {
    .capture_refuse(state, path, "It is a missing argument.",
                    "Give the argument a value in the function that makes the callback.")
  }
  value
}

.capture_check_size <- function(value, state, name) {
  size <- length(serialize(value, NULL))
  if (size <= .capture_limit()) {
    return(invisible())
  }
  biggest <- name
  if (length(state$bindings)) {
    sizes <- vapply(state$bindings, function(v) length(serialize(v, NULL)), numeric(1))
    biggest <- names(state$bindings)[[which.max(sizes)]]
    size <- max(sizes)
  }
  label <- state$label
  cli::cli_abort(
    c(
      "{label} captures {.var {biggest}}, which is {prettyNum(size, big.mark = ',')} bytes when serialised.",
      i = "Each job sends captured values to its worker. The limit is {prettyNum(.capture_limit(), big.mark = ',')} bytes.",
      i = "Read large inputs inside the function, or raise {.code options(bayesim.capture_limit = )}."
    ),
    call = NULL
  )
}

# The names that a closure reads from its environment, sorted.
.closure_names <- function(f) {
  globals <- suppressWarnings(codetools::findGlobals(f, merge = FALSE))
  names <- c(globals$functions, globals$variables, .free_reads(f))
  .radix_sort(unique(names))
}

# Every name that the code of `f` can read before it has surely assigned
# it. An assignment counts only where it surely runs before what follows:
# in a sequence, in both branches of an `if`, or as a loop variable.
.free_reads <- function(f) {
  acc <- new.env(parent = emptyenv())
  acc$reads <- character()
  read <- function(name, assigned) {
    if (nzchar(name) && !(name %in% assigned)) {
      acc$reads <- c(acc$reads, name)
    }
  }
  empty <- function(x) is.symbol(x) && !nzchar(as.character(x))
  each <- function(parts, assigned) {
    for (i in seq_along(parts)) {
      if (!empty(parts[[i]])) walk(parts[[i]], assigned)
    }
  }
  # The variable that a replacement such as `names(x$a) <- v` changes.
  target_var <- function(x) {
    while (is.call(x)) x <- x[[2L]]
    if (is.symbol(x) || is.character(x)) as.character(x)
  }
  formals_of <- function(fl, assigned) {
    inner <- union(assigned, names(fl))
    for (i in seq_along(fl)) {
      if (!empty(fl[[i]])) walk(fl[[i]], inner)
    }
    inner
  }
  walk <- function(x, assigned) {
    if (is.symbol(x)) {
      read(as.character(x), assigned)
      return(assigned)
    }
    if (!is.call(x)) {
      return(assigned)
    }
    head <- x[[1L]]
    fn <- if (is.symbol(head)) as.character(head) else ""
    parts <- as.list(x)[-1L]
    if (fn %in% c("<-", "=") && length(parts) == 2L) {
      assigned <- walk(parts[[2L]], assigned)
      target <- parts[[1L]]
      if (is.call(target)) {
        # A replacement reads the variable and the other arguments.
        walk(target, assigned)
      }
      var <- target_var(target)
      return(if (is.null(var)) assigned else union(assigned, var))
    }
    if (fn == "<<-" && length(parts) == 2L) {
      walk(parts[[2L]], assigned)
      walk(parts[[1L]], character())
      return(assigned)
    }
    if (fn %in% c("{", "(")) {
      for (i in seq_along(parts)) {
        if (!empty(parts[[i]])) assigned <- walk(parts[[i]], assigned)
      }
      return(assigned)
    }
    if (fn == "if") {
      assigned <- walk(parts[[1L]], assigned)
      yes <- walk(parts[[2L]], assigned)
      no <- if (length(parts) > 2L) walk(parts[[3L]], assigned) else assigned
      return(intersect(yes, no))
    }
    if (fn == "for") {
      var <- as.character(parts[[1L]])
      assigned <- walk(parts[[2L]], assigned)
      walk(parts[[3L]], union(assigned, var))
      return(union(assigned, var))
    }
    if (fn == "function") {
      inner <- formals_of(parts[[1L]], assigned)
      walk(parts[[2L]], inner)
      return(assigned)
    }
    if (fn == "quote") {
      return(assigned)
    }
    if (!is.symbol(head)) {
      walk(head, assigned)
    } else {
      read(fn, assigned)
    }
    # Arguments are promises: an assignment in one may never run.
    each(parts, assigned)
    assigned
  }
  assigned <- formals_of(formals(f), character())
  walk(body(f), assigned)
  unique(acc$reads)
}

# Whether the code of `f` uses `...` or `..1` where no function of its own
# declares `...`.
.uses_free_dots <- function(f) {
  dots <- function(x) {
    name <- as.character(x)
    name == "..." || grepl("^[.][.][0-9]+$", name)
  }
  empty <- function(x) is.symbol(x) && !nzchar(as.character(x))
  walk <- function(x, declared) {
    if (is.symbol(x)) {
      return(!declared && dots(x))
    }
    if (!is.call(x) && !is.pairlist(x)) {
      return(FALSE)
    }
    if (is.call(x) && identical(x[[1L]], as.name("function"))) {
      # A nested function declares its own dots; walk only its body.
      return(walk(x[[3L]], declared || "..." %in% names(x[[2L]])))
    }
    parts <- as.list(x)
    for (i in seq_along(parts)) {
      if (!empty(parts[[i]]) && walk(parts[[i]], declared)) {
        return(TRUE)
      }
    }
    FALSE
  }
  walk(body(f), "..." %in% names(formals(f)))
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

# Identity ---------------------------------------------------------------------
#
# The identity of a transport form is a digest of .stable(), a form of it
# that digests the same in every session:
#
# * Atomic values and lists carry their attributes in their order, classed
#   lists included. An S4 object carries its data part and its slots.
# * A closure counts by its formals and body, its environment, the packages
#   its body names with `::`, and its attributes. Its environment is a copy
#   that holds only what the closure reads, so a binding that it does not
#   read does not count.
# * A package function counts by package, version and name. Its body does
#   not count, because the version names it.
# * A primitive counts by name.
# * A formula counts by its text, its environment and its attributes.
# * An environment counts by its bindings, in radix order of their names,
#   its parent and its attributes. Namespaces and package environments count
#   by name and version.
# * Closures and environments that the walk has already entered count by
#   their place in the walk, so cycles end.
#
# Every order that feeds a digest is a radix order, which does not depend on
# the locale.

.value_id <- function(value) {
  .digest(list("value", .stable(value)))
}

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
    data <- x
    attributes(data) <- NULL
    data <- if (typeof(x) != "S4") .stable(data, state)
    return(list("S4", typeof(x), data, attrs))
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
  pkgs <- .radix_sort(unique(.qualified_packages(body(f))))
  list(
    "closure",
    .fn_shape(f),
    .stable_env(env, state),
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
    pkg <- .env_package(e)
    return(list("environment", .env_label(e), .package_version(pkg)))
  }
  seen <- .stable_enter(e, state)
  if (!is.null(seen)) {
    return(seen)
  }
  names <- .radix_sort(ls(e, all.names = TRUE, sorted = FALSE))
  values <- lapply(names, function(name) .stable(.get_binding(name, e), state))
  list(
    "environment",
    names,
    values,
    .stable_env(parent.env(e), state),
    isS4(e),
    .stable_attrs(attributes(e), state)
  )
}

.stable_formula <- function(x, state) {
  env <- environment(x)
  attrs <- attributes(x)
  attrs$.Environment <- NULL
  text <- x
  attributes(text) <- NULL
  list(
    "formula",
    deparse(text),
    if (is.environment(env)) .stable_env(env, state),
    .stable_attrs(attrs, state)
  )
}

# A binding's value. Capture refuses the bindings that get() cannot read;
# this keeps the identity defined for any environment.
.get_binding <- function(name, env) {
  if (bindingIsActive(name, env)) {
    return(list("active binding"))
  }
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
  for (name in .radix_sort(ls(env, all.names = TRUE, sorted = FALSE))) {
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

# Study capture ------------------------------------------------------------------

# Capture every callback, setting and condition of a study once, and compute
# the identities that keys are built from. The result is what workers
# receive.
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
    settings <- .capture(m$settings, label("Method \"", id, "\""), "settings")
    list(
      fit = fit$fn,
      extract = lapply(extract, `[[`, "fn"),
      settings = settings$value,
      fit_id = .digest(list(fit$id, settings$id, m$version)),
      extract_id = .digest(list(
        lapply(extract, `[[`, "id")[.radix_sort(names(extract))],
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
  conditions <- study$conditions
  rows <- lapply(seq_len(nrow(conditions)), function(i) {
    .capture(
      as.list(conditions[i, , drop = FALSE]),
      label("Condition \"", conditions$condition_id[[i]], "\"")
    )
  })
  list(
    name = study$name,
    generate = gen$fn,
    prepare = prep$fn,
    generate_id = .digest(list(gen$id, study$version)),
    prepare_id = .digest(list(prep$id, study$version)),
    methods = methods,
    measures = measures,
    comparisons = comparisons,
    conditions = lapply(rows, `[[`, "value"),
    condition_keys = vapply(rows, `[[`, "", "id"),
    retention = study$retention,
    compress = study$compress %||% "gzip",
    required = .sg_required(study)
  )
}
