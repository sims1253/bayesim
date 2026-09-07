# Experimental interface and reference executor. Source into an R session.
# This does not replace bayesim's supported simulation_config() engine.

study <- function(name, conditions, generate, prepare = NULL, version = "1") {
  .sg_name(name)
  if (
    !is.data.frame(conditions) ||
      !nrow(conditions) ||
      !is.character(conditions$condition_id) ||
      anyNA(conditions$condition_id) ||
      !all(nzchar(conditions$condition_id)) ||
      anyDuplicated(conditions$condition_id)
  ) {
    stop("conditions needs a unique, nonempty character condition_id column")
  }
  .sg_function(generate)
  if (!is.null(prepare)) {
    .sg_function(prepare)
  }
  .sg_name(version)
  structure(
    list(
      name = name,
      conditions = conditions,
      generate = generate,
      prepare = prepare,
      version = version,
      methods = list(),
      measures = list(),
      comparisons = list(),
      retention = character(),
      compress = "gzip"
    ),
    class = "bayesim_study_declaration"
  )
}

study_method <- function(
  fit,
  extract = list(),
  settings = list(),
  version = "1"
) {
  .sg_function(fit)
  .sg_named(extract)
  lapply(extract, .sg_function)
  if (any(names(extract) %in% c("fit", "data"))) {
    stop("fit and data are built-in artifacts; use other extractor names")
  }
  if (!is.list(settings)) {
    stop("settings must be a list")
  }
  .sg_named(settings)
  .sg_name(version)
  structure(
    list(fit = fit, extract = extract, settings = settings, version = version),
    class = "bayesim_study_method"
  )
}

with_methods <- function(study, ...) {
  .sg_study(study)
  methods <- list(...)
  .sg_named(methods)
  if (
    !length(methods) ||
      !all(vapply(methods, inherits, logical(1), "bayesim_study_method"))
  ) {
    stop("Supply named study_method() objects")
  }
  if (any(names(methods) %in% names(study$methods))) {
    stop("Method IDs already exist")
  }
  study$methods <- c(study$methods, methods)
  study
}

with_measure <- function(
  study,
  name,
  compute,
  needs = character(),
  version = "1"
) {
  .sg_study(study)
  .sg_name(name)
  .sg_function(compute)
  .sg_needs(needs)
  .sg_name(version)
  if (name %in% names(study$measures)) {
    stop("Measurement name already exists")
  }
  study$measures[[name]] <- list(
    compute = compute,
    needs = unique(needs),
    version = version
  )
  study
}

with_comparison <- function(
  study,
  name,
  compute,
  by = "dataset_id",
  needs = character(),
  version = "1"
) {
  .sg_study(study)
  .sg_name(name)
  .sg_function(compute)
  .sg_needs(needs)
  .sg_needs(by)
  .sg_name(version)
  if (name %in% names(study$comparisons)) {
    stop("Comparison name already exists")
  }
  if (length(needs) && !"dataset_id" %in% by) {
    stop("Artifact-consuming comparisons must group within dataset_id")
  }
  study$comparisons[[name]] <- list(
    compute = compute,
    by = by,
    needs = unique(needs),
    version = version
  )
  study
}

with_retention <- function(study, artifacts = character(), compress = "gzip") {
  .sg_study(study)
  .sg_needs(artifacts)
  if (!compress %in% c("gzip", "bzip2", "xz")) {
    stop("Use gzip, bzip2 or xz compression")
  }
  study$retention <- unique(artifacts)
  study$compress <- compress
  study
}

plan_study <- function(study, replicates, seed = 1L) {
  .sg_study(study)
  .sg_count(replicates, "replicates")
  .sg_count(seed, "seed", zero = TRUE)
  if (!length(study$methods)) {
    stop("Add at least one method")
  }
  needs <- .sg_required(study)
  by <- unique(unlist(lapply(study$comparisons, `[[`, "by"), use.names = FALSE))
  columns <- c(
    "condition_id",
    "dataset_id",
    "replicate",
    "method_id",
    paste0("condition_", setdiff(names(study$conditions), "condition_id")),
    paste0(
      "method_",
      unique(unlist(lapply(study$methods, function(x) names(x$settings))))
    )
  )
  if (length(setdiff(by, columns))) {
    stop(
      "Unknown comparison grouping columns: ",
      paste(setdiff(by, columns), collapse = ", ")
    )
  }
  for (id in names(study$methods)) {
    available <- c("data", "fit", names(study$methods[[id]]$extract))
    missing <- setdiff(needs, available)
    if (length(missing)) {
      stop("Method '", id, "' cannot supply: ", paste(missing, collapse = ", "))
    }
  }
  dependencies <- .sg_bind(lapply(names(study$measures), function(id) {
    data.frame(
      stage = "per fit",
      name = id,
      needs = paste(study$measures[[id]]$needs, collapse = ", "),
      stringsAsFactors = FALSE
    )
  }))
  dependencies <- .sg_bind(c(
    list(dependencies),
    lapply(names(study$comparisons), function(id) {
      x <- study$comparisons[[id]]
      data.frame(
        stage = if (length(x$needs)) {
          "within dataset, before eviction"
        } else {
          "saved measurements"
        },
        name = id,
        needs = paste(x$needs, collapse = ", "),
        stringsAsFactors = FALSE
      )
    })
  ))
  structure(
    list(
      study = study$name,
      conditions = nrow(study$conditions),
      replicates = as.integer(replicates),
      datasets = nrow(study$conditions) * as.double(replicates),
      methods = names(study$methods),
      fits = nrow(study$conditions) *
        as.double(replicates) *
        length(study$methods),
      prepare = !is.null(study$prepare),
      required = needs,
      retained = study$retention,
      discarded = setdiff(needs, study$retention),
      dependencies = dependencies,
      seed = as.integer(seed)
    ),
    class = "bayesim_study_plan"
  )
}

print.bayesim_study_plan <- function(x, ...) {
  cat(
    x$study,
    "\n",
    x$conditions,
    " conditions x ",
    x$replicates,
    " replicates = ",
    format(x$datasets, scientific = FALSE),
    " datasets\n",
    length(x$methods),
    " methods; ",
    format(x$fits, scientific = FALSE),
    " fits\n",
    sep = ""
  )
  if (x$prepare) {
    cat("Prepare once per condition.\n")
  }
  print(x$dependencies, row.names = FALSE)
  cat(
    "Retain: ",
    if (length(x$retained)) {
      paste(x$retained, collapse = ", ")
    } else {
      "measurements only"
    },
    "\n",
    sep = ""
  )
  cat(
    "Discard after use: ",
    paste(x$discarded, collapse = ", "),
    "\n",
    sep = ""
  )
  cat("Plan only: no data, models or cache records have been evaluated.\n")
  invisible(x)
}

# Dependencies are function bodies, formals, explicit versions and settings.
# External files, globals and package upgrades require a version change.
.sg_hash <- function(x) digest::digest(x, algo = "sha256")
.sg_recipe <- function(fn, version) {
  if (is.null(fn)) {
    return(NULL)
  }
  list(formals = formals(fn), body = body(fn), version = version)
}
.sg_name <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("Expected one nonempty name or version")
  }
}
.sg_function <- function(x) if (!is.function(x)) stop("Expected a function")
.sg_needs <- function(x) {
  if (!is.character(x) || anyNA(x) || !all(nzchar(x)) || anyDuplicated(x)) {
    stop("Expected unique nonempty character names")
  }
}
.sg_named <- function(x) {
  if (
    length(x) &&
      (is.null(names(x)) ||
        anyNA(names(x)) ||
        !all(nzchar(names(x))) ||
        anyDuplicated(names(x)))
  ) {
    stop("Every entry must have a unique nonempty name")
  }
}
.sg_study <- function(x) {
  if (!inherits(x, "bayesim_study_declaration")) {
    stop("Expected a study() declaration")
  }
}
.sg_count <- function(x, name, zero = FALSE) {
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      !is.finite(x) ||
      x != floor(x) ||
      x < (if (zero) 0 else 1) ||
      x > .Machine$integer.max
  ) {
    stop(name, " must be a whole number in the supported integer range")
  }
}
.sg_required <- function(study) {
  unique(c(
    study$retention,
    unlist(
      lapply(c(study$measures, study$comparisons), `[[`, "needs"),
      use.names = FALSE
    )
  ))
}
.sg_bind <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, rows)
  if (!length(rows)) {
    return(data.frame())
  }
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  list_columns <- columns[vapply(
    columns,
    function(nm) {
      any(vapply(
        rows,
        function(x) nm %in% names(x) && is.list(x[[nm]]),
        logical(1)
      ))
    },
    logical(1)
  )]
  rows <- lapply(rows, function(x) {
    for (nm in setdiff(columns, names(x))) {
      x[[nm]] <- if (nm %in% list_columns) {
        rep(list(NULL), nrow(x))
      } else {
        rep(NA, nrow(x))
      }
    }
    x[, columns, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
.sg_words <- function(key) {
  h <- .sg_hash(key)
  vapply(
    seq_len(6L),
    function(i) strtoi(substr(h, (i - 1L) * 7L + 1L, i * 7L), 16L) + 1L,
    integer(1)
  )
}
.sg_invoke <- function(fn, args, key) {
  withr::local_seed(1L, .rng_kind = "L'Ecuyer-CMRG")
  stream <- get(".Random.seed", envir = .GlobalEnv)
  stream[2:7] <- .sg_words(key)
  assign(".Random.seed", stream, envir = .GlobalEnv)
  if ("context" %in% names(args)) {
    args$context$seed <- .sg_words(key)[1L]
  }
  do.call(fn, args)
}
.sg_path <- function(path, kind, key) {
  if (is.null(path)) {
    return(NULL)
  }
  file.path(path, kind, substr(key, 1, 2), paste0(key, ".rds"))
}
.sg_read <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }
  record <- tryCatch(readRDS(path), error = function(e) {
    stop("Unreadable experimental cache record: ", path)
  })
  if (!is.list(record) || !identical(record$checksum, .sg_hash(record$value))) {
    stop("Checksum mismatch: ", path)
  }
  record$value
}
.sg_write <- function(value, path, compress) {
  if (is.null(path)) {
    return(invisible(NULL))
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".pending-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(
    list(value = value, checksum = .sg_hash(value)),
    tmp,
    compress = compress
  )
  if (!file.rename(tmp, path)) {
    stop("Could not commit experimental record: ", path)
  }
  invisible(NULL)
}
.sg_open <- function(study, seed, path) {
  if (is.null(path)) {
    return(NULL)
  }
  if (
    !is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)
  ) {
    stop("path must be NULL or a directory")
  }
  manifest <- file.path(path, "experimental-study.rds")
  expected <- list(format = 1L, study = study$name, seed = as.integer(seed))
  if (file.exists(manifest)) {
    if (!identical(.sg_read(manifest), expected)) {
      stop("Study name, seed or experimental format differs at this path")
    }
  } else {
    if (file.exists(path) && !dir.exists(path)) {
      stop("path is a file")
    }
    if (
      dir.exists(path) &&
        length(list.files(path, all.files = TRUE, no.. = TRUE))
    ) {
      stop("path contains unrelated files")
    }
    .sg_write(expected, manifest, study$compress)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
.sg_condition <- function(study, id) {
  hit <- match(id, study$conditions$condition_id)
  if (is.na(hit)) {
    stop("Unknown condition_id: ", id)
  }
  as.list(study$conditions[hit, , drop = FALSE])
}
.sg_prepare <- function(study, condition, seed, path) {
  key <- .sg_hash(list(
    study$name,
    seed,
    condition,
    .sg_recipe(study$prepare, study$version)
  ))
  context <- list(study = study$name, condition = condition, seed = seed)
  if (is.null(study$prepare)) {
    return(list(value = NULL, key = key))
  }
  file <- .sg_path(path, "prepared", key)
  value <- .sg_read(file)
  if (is.null(value)) {
    value <- list(
      value = .sg_invoke(
        study$prepare,
        list(condition = condition, context = context),
        key
      ),
      key = key
    )
    .sg_write(value, file, study$compress)
  }
  value
}
.sg_generate <- function(study, condition, context, key) {
  context$method_id <- NULL
  context$settings <- NULL
  .sg_invoke(
    study$generate,
    list(condition = condition, context = context),
    key
  )
}
.sg_metadata <- function(context) {
  x <- list(
    condition_id = context$condition$condition_id,
    dataset_id = context$dataset_id,
    replicate = context$replicate,
    method_id = context$method_id
  )
  for (prefix in c("condition", "method")) {
    values <- if (prefix == "condition") context$condition else context$settings
    for (nm in setdiff(names(values), "condition_id")) {
      if (is.atomic(values[[nm]]) && length(values[[nm]]) == 1L) {
        x[[paste0(prefix, "_", nm)]] <- values[[nm]]
      }
    }
  }
  as.data.frame(x, stringsAsFactors = FALSE)
}
.sg_rows <- function(rows, context, name) {
  if (!is.data.frame(rows)) {
    stop("Measurement '", name, "' must return a data frame")
  }
  if (!nrow(rows)) {
    return(data.frame())
  }
  meta <- .sg_metadata(context)
  if (any(names(rows) %in% c(names(meta), "measurement"))) {
    stop("Measurement output uses reserved identity columns")
  }
  cbind(
    meta[rep(1L, nrow(rows)), , drop = FALSE],
    measurement = name,
    rows,
    row.names = NULL
  )
}
.sg_artifacts <- function(method, fit, data, context, needs, key) {
  out <- list()
  for (name in needs) {
    out[name] <- list(
      if (name == "fit") {
        fit
      } else if (name == "data") {
        data
      } else {
        .sg_invoke(
          method$extract[[name]],
          list(fit = fit, data = data, context = context),
          list(key, "extract", name)
        )
      }
    )
  }
  out
}
.sg_measure <- function(study, artifacts, context, cached, key) {
  rows <- list()
  for (name in names(study$measures)) {
    spec <- study$measures[[name]]
    fingerprint <- .sg_hash(list(
      .sg_recipe(spec$compute, spec$version),
      spec$needs
    ))
    existing <- cached[[name]]
    if (!is.null(existing) && identical(existing$fingerprint, fingerprint)) {
      rows[[name]] <- existing
    } else {
      missing <- setdiff(spec$needs, names(artifacts))
      if (length(missing)) {
        stop(
          "Measurement '",
          name,
          "' needs discarded artifacts: ",
          paste(missing, collapse = ", "),
          ". Use a new run path to refit with the required retention."
        )
      }
      result <- .sg_invoke(
        spec$compute,
        list(artifacts = artifacts[spec$needs], context = context),
        list(key, "measure", name, fingerprint)
      )
      rows[[name]] <- list(
        fingerprint = fingerprint,
        rows = .sg_rows(result, context, name)
      )
    }
  }
  rows
}
.sg_compare <- function(
  study,
  results,
  attempts,
  bundles = list(),
  raw = FALSE,
  path = NULL
) {
  output <- list()
  for (name in names(study$comparisons)) {
    spec <- study$comparisons[[name]]
    if ((length(spec$needs) > 0L) != raw) {
      next
    }
    # Group from attempts so entirely failed groups remain observable.
    grouping <- attempts[, spec$by, drop = FALSE]
    groups <- if (!length(spec$by)) {
      rep(1L, nrow(attempts))
    } else {
      getFromNamespace("group_ids", "bayesim")(grouping, spec$by)
    }
    for (group in unique(groups)) {
      members <- attempts[groups == group, , drop = FALSE]
      selected <- if (nrow(results)) {
        results$dataset_id %in%
          members$dataset_id &
          results$method_id %in% members$method_id
      } else {
        logical()
      }
      # Match full dataset/method pairs; separate membership sets can cross pairs.
      if (nrow(results)) {
        selected <- paste(results$dataset_id, results$method_id, sep = "/") %in%
          paste(members$dataset_id, members$method_id, sep = "/")
      }
      rows <- results[selected, , drop = FALSE]
      context <- list(
        study = study$name,
        group = if (length(spec$by)) {
          as.list(grouping[which(groups == group)[1L], , drop = FALSE])
        } else {
          list()
        },
        attempts = members
      )
      key <- list(
        "comparison",
        name,
        .sg_recipe(spec$compute, spec$version),
        spec$needs,
        context$group,
        members[, setdiff(names(members), "elapsed"), drop = FALSE],
        rows
      )
      file <- if (raw) .sg_path(path, "comparisons", .sg_hash(key)) else NULL
      cached <- .sg_read(file)
      artifacts <- bundles[intersect(members$method_id, names(bundles))]
      if (raw && is.null(cached)) {
        for (id in members$method_id[members$status == "success"]) {
          if (length(setdiff(spec$needs, names(artifacts[[id]])))) {
            stop(
              "Comparison '",
              name,
              "' needs discarded artifacts for method '",
              id,
              "'. Refit in a new path with those artifacts retained."
            )
          }
        }
        artifacts <- lapply(artifacts, function(x) x[spec$needs])
      }
      value <- if (!is.null(cached)) {
        cached$value
      } else {
        .sg_invoke(
          spec$compute,
          list(results = rows, artifacts = artifacts, context = context),
          key
        )
      }
      if (!is.data.frame(value)) {
        stop("Comparison must return a data frame: ", name)
      }
      if (raw && is.null(cached)) {
        .sg_write(list(value = value), file, study$compress)
      }
      if (nrow(value)) {
        if (any(names(value) %in% c("comparison", spec$by))) {
          stop("Comparison output uses reserved grouping columns")
        }
        meta <- if (length(spec$by)) {
          grouping[rep(which(groups == group)[1L], nrow(value)), , drop = FALSE]
        } else {
          data.frame(row.names = seq_len(nrow(value)))
        }
        output[[length(output) + 1L]] <- cbind(
          meta,
          comparison = name,
          value,
          row.names = NULL
        )
      }
    }
  }
  .sg_bind(output)
}

# Dataset-level evaluator. Durable outcomes are committed after each fit.
evaluate_replicate <- function(
  study,
  condition_id,
  replicate,
  seed = 1L,
  path = NULL,
  .prepared = NULL
) {
  plan_study(study, replicate, seed)
  path <- .sg_open(study, seed, path)
  condition <- .sg_condition(study, condition_id)
  prepared <- if (is.null(.prepared)) {
    .sg_prepare(study, condition, seed, path)
  } else {
    .prepared
  }
  data_key <- .sg_hash(list(
    study$name,
    seed,
    condition,
    replicate,
    prepared$key,
    .sg_recipe(study$generate, study$version)
  ))
  context <- list(
    study = study$name,
    condition = condition,
    dataset_id = data_key,
    replicate = as.integer(replicate),
    prepared = prepared$value,
    method_id = "",
    settings = list(),
    seed = seed
  )
  data_file <- .sg_path(path, "datasets", data_key)
  saved_data <- .sg_read(data_file)
  data <- if (is.null(saved_data)) NULL else saved_data$data
  has_data <- !is.null(saved_data)
  needs <- .sg_required(study)
  attempts <- measurements <- bundles <- list()
  for (id in names(study$methods)) {
    method <- study$methods[[id]]
    context$method_id <- id
    context$settings <- method$settings
    fit_key <- .sg_hash(list(
      data_key,
      id,
      method$settings,
      .sg_recipe(method$fit, method$version)
    ))
    file <- .sg_path(path, "fits", fit_key)
    cached <- .sg_read(file)
    meta <- .sg_metadata(context)
    meta$fit_id <- fit_key
    if (!is.null(cached)) {
      if (
        !identical(
          cached$extract,
          lapply(method$extract, .sg_recipe, version = method$version)
        )
      ) {
        stop(
          "Extractor recipe changed for cached method '",
          id,
          "'. Change the method version to refit."
        )
      }
      if (cached$attempt$status != "success") {
        attempts[[id]] <- cached$attempt
        next
      }
      artifacts <- cached$artifacts
      if ("data" %in% needs && !"data" %in% names(artifacts)) {
        if (!has_data) {
          data <- .sg_generate(study, condition, context, data_key)
          has_data <- TRUE
        }
        artifacts["data"] <- list(data)
      }
      missing_retention <- setdiff(study$retention, names(artifacts))
      if (length(missing_retention)) {
        stop(
          "Requested retention needs discarded artifacts: ",
          paste(missing_retention, collapse = ", ")
        )
      }
      rows <- .sg_measure(
        study,
        artifacts,
        context,
        cached$measurements,
        fit_key
      )
      updated <- cached
      updated$measurements <- rows
      if (!identical(rows, cached$measurements)) {
        .sg_write(updated, file, study$compress)
      }
      attempts[[id]] <- cached$attempt
    } else {
      if (!has_data) {
        data <- .sg_generate(study, condition, context, data_key)
        has_data <- TRUE
        if ("data" %in% study$retention) {
          .sg_write(list(data = data), data_file, study$compress)
        }
      }
      started <- proc.time()[["elapsed"]]
      fit <- tryCatch(
        list(
          value = .sg_invoke(
            method$fit,
            list(data = data, context = context),
            list(fit_key, "fit")
          )
        ),
        error = function(e) list(error = conditionMessage(e))
      )
      attempt <- cbind(
        meta,
        status = if (is.null(fit$error)) "success" else "failed",
        error = if (is.null(fit$error)) NA_character_ else fit$error,
        elapsed = proc.time()[["elapsed"]] - started,
        stringsAsFactors = FALSE
      )
      if (!is.null(fit$error)) {
        .sg_write(
          list(
            attempt = attempt,
            artifacts = list(),
            measurements = list(),
            extract = lapply(
              method$extract,
              .sg_recipe,
              version = method$version
            )
          ),
          file,
          study$compress
        )
        attempts[[id]] <- attempt
        next
      }
      artifacts <- .sg_artifacts(
        method,
        fit$value,
        data,
        context,
        needs,
        fit_key
      )
      rows <- .sg_measure(study, artifacts, context, list(), fit_key)
      .sg_write(
        list(
          attempt = attempt,
          artifacts = artifacts[intersect(
            setdiff(
              unique(c(
                study$retention,
                unlist(
                  lapply(study$comparisons, `[[`, "needs"),
                  use.names = FALSE
                )
              )),
              "data"
            ),
            names(artifacts)
          )],
          measurements = rows,
          extract = lapply(method$extract, .sg_recipe, version = method$version)
        ),
        file,
        study$compress
      )
      attempts[[id]] <- attempt
      rm(fit)
    }
    measurements[[id]] <- .sg_bind(lapply(rows, `[[`, "rows"))
    group_needs <- unique(unlist(
      lapply(study$comparisons, `[[`, "needs"),
      use.names = FALSE
    ))
    bundles[[id]] <- artifacts[intersect(group_needs, names(artifacts))]
    rm(artifacts)
  }
  attempts <- .sg_bind(attempts)
  measurements <- .sg_bind(measurements)
  # Group artifacts live only until all methods for this dataset are available.
  comparisons <- .sg_compare(
    study,
    measurements,
    attempts,
    bundles,
    raw = TRUE,
    path = path
  )
  if (
    !is.null(path) && length(unlist(lapply(study$comparisons, `[[`, "needs")))
  ) {
    for (key in attempts$fit_id[attempts$status == "success"]) {
      file <- .sg_path(path, "fits", key)
      record <- .sg_read(file)
      extra <- setdiff(names(record$artifacts), study$retention)
      if (length(extra)) {
        record$artifacts[extra] <- NULL
        .sg_write(record, file, study$compress)
      }
    }
  }
  list(
    attempts = attempts,
    measurements = measurements,
    comparisons = comparisons
  )
}

run_study <- function(study, replicates, seed = 1L, path = NULL, workers = 1L) {
  plan <- plan_study(study, replicates, seed)
  .sg_count(workers, "workers")
  path <- .sg_open(study, seed, path)
  prepared <- lapply(study$conditions$condition_id, function(id) {
    .sg_prepare(study, .sg_condition(study, id), seed, path)
  })
  names(prepared) <- study$conditions$condition_id
  # Batches bound active artifacts; this reference runner still collects all scalar results.
  jobs <- expand.grid(
    condition_id = study$conditions$condition_id,
    replicate = seq_len(replicates),
    stringsAsFactors = FALSE
  )
  output <- vector("list", nrow(jobs))
  cl <- NULL
  if (workers > 1L) {
    cl <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    env <- environment(run_study)
    parallel::clusterExport(
      cl,
      c(
        ls(env, pattern = "^\\.sg_", all.names = TRUE),
        "evaluate_replicate",
        "plan_study"
      ),
      envir = env
    )
  }
  worker <- function(job, study, seed, path, prepared) {
    evaluate_replicate(
      study,
      job$condition_id,
      job$replicate,
      seed,
      path,
      .prepared = prepared[[job$condition_id]]
    )
  }
  for (start in seq.int(1L, nrow(jobs), by = workers)) {
    indices <- seq.int(start, min(nrow(jobs), start + workers - 1L))
    batch <- lapply(indices, function(i) as.list(jobs[i, , drop = FALSE]))
    values <- if (is.null(cl)) {
      lapply(batch, worker, study, seed, path, prepared)
    } else {
      parallel::parLapplyLB(cl, batch, worker, study, seed, path, prepared)
    }
    output[indices] <- values
  }
  run <- structure(
    list(
      study = study,
      plan = plan,
      path = path,
      attempts = .sg_bind(lapply(output, `[[`, "attempts")),
      measurements = .sg_bind(lapply(output, `[[`, "measurements")),
      comparisons = .sg_bind(lapply(output, `[[`, "comparisons"))
    ),
    class = "bayesim_study_run"
  )
  run$comparisons <- .sg_bind(list(
    run$comparisons,
    .sg_compare(study, run$measurements, run$attempts)
  ))
  run
}

assess_study <- function(run, policy = NULL) {
  if (!inherits(run, "bayesim_study_run")) {
    stop("Expected a run_study() result")
  }
  rows <- run$measurements
  keep <- if (is.null(policy)) {
    rep(TRUE, nrow(rows))
  } else {
    policy(rows, run$attempts)
  }
  if (!is.logical(keep) || length(keep) != nrow(rows) || anyNA(keep)) {
    stop("policy must return one nonmissing logical per measurement row")
  }
  selected <- rows[keep, , drop = FALSE]
  list(
    measurements = selected,
    excluded = rows[!keep, , drop = FALSE],
    attempts = run$attempts,
    comparisons = .sg_compare(run$study, selected, run$attempts),
    execution_comparisons = run$comparisons[
      if (nrow(run$comparisons)) {
        run$comparisons$comparison %in%
          names(Filter(function(x) length(x$needs) > 0L, run$study$comparisons))
      } else {
        logical()
      },
      ,
      drop = FALSE
    ]
  )
}
