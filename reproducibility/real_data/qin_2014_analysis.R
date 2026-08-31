#!/usr/bin/env Rscript

###############################################################################
# Qin et al. (2014) real-data analysis
#
# Portable entrypoint for the QinN_2014 liver-cirrhosis microbiome data in
# curatedMetagenomicData.  This script is deliberately source-safe: sourcing it
# defines functions but does not download data or run an analysis.
#
# Negative controls use the shared simulation implementation:
#   * response-free two-fold cross-fitting;
#   * a Gaussian regression on the first five raw-abundance PCs of X_{-j};
#   * stochastic conditional draws clipped to the abundance support; and
#   * direct replacement of feature j, with no ILR/CLR transformation and no
#     closure adjustment of the remaining features.
###############################################################################

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

qin_script_path <- function() {
  # When this file is sourced, an `ofile` frame identifies this file even if
  # commandArgs() refers to a different parent script.
  frame_count <- sys.nframe()
  if (frame_count > 0L) {
    for (i in rev(seq_len(frame_count))) {
      candidate <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
      if (!is.null(candidate) && nzchar(candidate)) {
        return(normalizePath(candidate, mustWork = FALSE))
      }
    }
  }

  # Unlike source(), sys.source() does not populate an `ofile` frame. Recover
  # its path from the active call so this file remains safe to load into a
  # dedicated environment during tests and downstream reuse.
  calls <- sys.calls()
  for (i in rev(seq_along(calls))) {
    call_i <- calls[[i]]
    if (length(call_i) < 2L || !identical(as.character(call_i[[1L]]), "sys.source")) {
      next
    }
    candidate <- tryCatch(
      eval(call_i[[2L]], envir = sys.frame(max(i - 1L, 0L))),
      error = function(e) NULL
    )
    if (is.character(candidate) && length(candidate) == 1L && file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE))
  }

  candidate <- file.path(
    getwd(), "reproducibility", "real_data", "qin_2014_analysis.R"
  )
  if (file.exists(candidate)) {
    return(normalizePath(candidate, mustWork = TRUE))
  }

  stop(
    "Cannot determine the location of qin_2014_analysis.R. ",
    "Run it with Rscript or source the file directly."
  )
}

QIN_SCRIPT_FILE <- qin_script_path()
QIN_REPO_ROOT <- normalizePath(
  file.path(dirname(QIN_SCRIPT_FILE), "..", ".."),
  mustWork = FALSE
)
QIN_SHARED_SCRIPT <- file.path(
  QIN_REPO_ROOT,
  "reproducibility",
  "simulation",
  "run_simulations.R"
)
QIN_WORKFLOW_VERSION <- "2026-08-31-qin-v1"

qin_is_direct_execution <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) {
    return(FALSE)
  }

  command_file <- normalizePath(
    sub("^--file=", "", file_arg),
    mustWork = FALSE
  )
  identical(command_file, normalizePath(QIN_SCRIPT_FILE, mustWork = FALSE))
}

qin_usage <- function() {
  paste(
    "Usage:",
    "  Rscript reproducibility/real_data/qin_2014_analysis.R [options]",
    "",
    "Options:",
    "  --output PATH       Output directory (default:",
    "                      reproducibility/results/qin_2014)",
    "  --q VALUES          Comma-separated FDR levels (default: 0.05,0.10)",
    "  --splits N          Number of MDS splits (default: 50)",
    "  --cores N           Threads passed to XGBoost (default: 1)",
    "  --seed N            Master seed (default: 20260824)",
    "  --pc N              Number of raw-abundance PCs (default: 5)",
    "  --prevalence X      Minimum taxon prevalence (default: 0.05)",
    "  --dry-run           Validate configuration without downloading data",
    "  --verbose           Print method progress",
    "  --help              Show this message",
    "",
    "Environment overrides:",
    "  FDR_QIN_OUTPUT, FDR_QIN_Q, FDR_QIN_SPLITS, FDR_QIN_CORES,",
    "  FDR_QIN_SEED, FDR_QIN_PC, FDR_QIN_PREVALENCE",
    sep = "\n"
  )
}

parse_flag_value <- function(args, i, option) {
  token <- args[[i]]
  prefix <- paste0(option, "=")

  if (startsWith(token, prefix)) {
    return(list(value = substring(token, nchar(prefix) + 1L), consumed = 1L))
  }

  if (identical(token, option)) {
    if (i == length(args)) {
      stop("Missing value after ", option, ".")
    }
    return(list(value = args[[i + 1L]], consumed = 2L))
  }

  NULL
}

parse_q_values <- function(x) {
  values <- suppressWarnings(as.numeric(strsplit(x, ",", fixed = TRUE)[[1L]]))
  if (length(values) < 1L || any(!is.finite(values)) ||
      any(values <= 0 | values >= 1)) {
    stop("--q must contain comma-separated values strictly between 0 and 1.")
  }
  sort(unique(values))
}

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

parse_qin_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- list(
    output = env_or_default(
      "FDR_QIN_OUTPUT",
      file.path("reproducibility", "results", "qin_2014")
    ),
    q_values = parse_q_values(env_or_default("FDR_QIN_Q", "0.05,0.10")),
    n_splits = as.numeric(env_or_default("FDR_QIN_SPLITS", "50")),
    cores = as.numeric(env_or_default("FDR_QIN_CORES", "1")),
    seed = as.numeric(env_or_default("FDR_QIN_SEED", "20260824")),
    K = as.numeric(env_or_default("FDR_QIN_PC", "5")),
    prevalence = as.numeric(env_or_default("FDR_QIN_PREVALENCE", "0.05")),
    dry_run = FALSE,
    verbose = FALSE,
    help = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]

    if (token %in% c("--help", "-h")) {
      config$help <- TRUE
      i <- i + 1L
      next
    }
    if (identical(token, "--dry-run")) {
      config$dry_run <- TRUE
      i <- i + 1L
      next
    }
    if (identical(token, "--verbose")) {
      config$verbose <- TRUE
      i <- i + 1L
      next
    }

    matched <- FALSE
    specifications <- list(
      output = "--output",
      q_values = "--q",
      n_splits = "--splits",
      cores = "--cores",
      seed = "--seed",
      K = "--pc",
      prevalence = "--prevalence"
    )

    for (field in names(specifications)) {
      parsed <- parse_flag_value(args, i, specifications[[field]])
      if (!is.null(parsed)) {
        value <- parsed$value
        if (field == "q_values") {
          value <- parse_q_values(value)
        } else if (field %in% c("n_splits", "cores", "seed", "K")) {
          value <- suppressWarnings(as.numeric(value))
        } else if (field == "prevalence") {
          value <- suppressWarnings(as.numeric(value))
        }
        config[[field]] <- value
        i <- i + parsed$consumed
        matched <- TRUE
        break
      }
    }

    if (!matched) {
      stop("Unknown option: ", token, "\n\n", qin_usage())
    }
  }

  if (!is.character(config$output) || length(config$output) != 1L ||
      !nzchar(config$output)) {
    stop("Output directory must be a non-empty path.")
  }
  integer_fields <- c("n_splits", "cores", "seed", "K")
  invalid_integer <- vapply(
    config[integer_fields],
    function(value) {
      length(value) != 1L || !is.finite(value) || value != floor(value) ||
        value < 1 || value > .Machine$integer.max
    },
    logical(1L)
  )
  if (any(invalid_integer)) {
    stop(
      "These options must be positive integers: ",
      paste(integer_fields[invalid_integer], collapse = ", "),
      "."
    )
  }
  config[integer_fields] <- lapply(config[integer_fields], as.integer)
  if (!is.finite(config$n_splits) || config$n_splits < 1L) {
    stop("--splits must be a positive integer.")
  }
  if (!is.finite(config$cores) || config$cores < 1L) {
    stop("--cores must be a positive integer.")
  }
  if (!is.finite(config$seed)) {
    stop("--seed must be an integer.")
  }
  if (!is.finite(config$K) || config$K < 1L) {
    stop("--pc must be a positive integer.")
  }
  if (!is.finite(config$prevalence) ||
      config$prevalence < 0 || config$prevalence > 1) {
    stop("--prevalence must be between 0 and 1.")
  }

  config$output <- path.expand(config$output)
  if (!grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", config$output)) {
    config$output <- file.path(QIN_REPO_ROOT, config$output)
  }
  config$output <- normalizePath(config$output, mustWork = FALSE)
  config
}

check_qin_packages <- function() {
  required <- c(
    "curatedMetagenomicData",
    "SummarizedExperiment",
    "glmnet",
    "xgboost",
    "knockoff"
  )
  missing <- required[
    !vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1L))
  ]

  if (length(missing) > 0L) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "), ".\n",
      "Install dependencies as described in README.md; this script never ",
      "installs packages automatically."
    )
  }
  invisible(required)
}

load_shared_simulation_api <- function() {
  if (!file.exists(QIN_SHARED_SCRIPT)) {
    stop("Shared simulation code was not found at: ", QIN_SHARED_SCRIPT)
  }

  api <- new.env(parent = globalenv())
  sys.source(QIN_SHARED_SCRIPT, envir = api)

  required_functions <- c(
    "generate_gaussian_pc_controls",
    "run_delporte_mds",
    "run_dai_mds",
    "run_univariate_by",
    "run_knockoffs",
    "default_xgboost_control"
  )
  missing <- required_functions[
    !vapply(
      required_functions,
      exists,
      logical(1L),
      envir = api,
      mode = "function",
      inherits = FALSE
    )
  ]

  if (length(missing) > 0L) {
    stop(
      "The shared simulation API is missing: ",
      paste(missing, collapse = ", "), "."
    )
  }
  api
}

as_numeric_measure <- function(x) {
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  text <- trimws(as.character(x))
  text[text %in% c("", "NA", "NaN", "N/A")] <- NA_character_
  suppressWarnings(as.numeric(sub(
    "^.*?(-?[0-9]+(?:\\.[0-9]+)?).*$",
    "\\1",
    text,
    perl = TRUE
  )))
}

extract_packed_biomarker <- function(names_column, values_column, pattern) {
  out <- rep(NA_real_, length(names_column))
  for (i in seq_along(names_column)) {
    if (is.na(names_column[[i]]) || is.na(values_column[[i]])) {
      next
    }
    biomarker_names <- trimws(strsplit(
      as.character(names_column[[i]]), ";", fixed = TRUE
    )[[1L]])
    biomarker_values <- trimws(strsplit(
      as.character(values_column[[i]]), ";", fixed = TRUE
    )[[1L]])
    hit <- grep(pattern, biomarker_names, ignore.case = TRUE)
    if (length(hit) > 0L && hit[[1L]] <= length(biomarker_values)) {
      out[[i]] <- as_numeric_measure(biomarker_values[[hit[[1L]]]])
    }
  }
  out
}

extract_albumin <- function(metadata) {
  direct_names <- names(metadata)[
    tolower(names(metadata)) %in% c("albumin", "albumine")
  ]
  if (length(direct_names) > 0L) {
    return(as_numeric_measure(metadata[[direct_names[[1L]]]]))
  }

  approximate <- grep("albumin", names(metadata), ignore.case = TRUE, value = TRUE)
  if (length(approximate) > 0L) {
    return(as_numeric_measure(metadata[[approximate[[1L]]]]))
  }

  if (all(c("biomarker_name", "biomarker_value") %in% names(metadata))) {
    return(extract_packed_biomarker(
      metadata$biomarker_name,
      metadata$biomarker_value,
      pattern = "albumin"
    ))
  }

  stop(
    "Could not identify albumin in QinN_2014 metadata. Available columns: ",
    paste(names(metadata), collapse = ", ")
  )
}

extract_cirrhosis <- function(metadata) {
  outcome <- rep(NA_real_, nrow(metadata))

  if ("control" %in% names(metadata)) {
    status <- tolower(trimws(as.character(metadata$control)))
    outcome[status == "case"] <- 1
    outcome[grepl("control|healthy", status)] <- 0
  } else if ("study_condition" %in% names(metadata)) {
    status <- tolower(trimws(as.character(metadata$study_condition)))
    outcome[grepl("cirrhosis", status)] <- 1
    outcome[status %in% c("control", "healthy")] <- 0
  } else if ("disease" %in% names(metadata)) {
    status <- tolower(trimws(as.character(metadata$disease)))
    outcome[grepl("cirrhosis", status)] <- 1
    outcome[status %in% c("control", "healthy")] <- 0
  } else {
    stop(
      "Could not identify cirrhosis/control status. Available columns: ",
      paste(names(metadata), collapse = ", ")
    )
  }

  outcome
}

qin_object_md5 <- function(object) {
  path <- tempfile("qin_input_", fileext = ".rds")
  on.exit(if (file.exists(path)) unlink(path), add = TRUE)
  saveRDS(object, path, version = 3L, compress = FALSE)
  unname(tools::md5sum(path))
}

load_qin_2014 <- function(prevalence_threshold = 0.05) {
  message("Downloading/loading QinN_2014.relative_abundance ...")
  experiments <- curatedMetagenomicData::curatedMetagenomicData(
    # ExperimentHub resource titles carry a release-date prefix.
    "QinN_2014[.]relative_abundance$",
    dryrun = FALSE,
    counts = FALSE,
    rownames = "short"
  )
  if (length(experiments) < 1L) {
    stop("QinN_2014.relative_abundance was not returned.")
  }
  if (length(experiments) > 1L) {
    warning("More than one QinN_2014 resource was returned; using the first.")
  }

  resource_names <- names(experiments)
  selected_resource <- if (length(resource_names) >= 1L &&
                           !is.na(resource_names[[1L]]) &&
                           nzchar(resource_names[[1L]])) {
    resource_names[[1L]]
  } else {
    "unnamed_first_resource"
  }
  experiment <- experiments[[1L]]
  assay <- as.matrix(SummarizedExperiment::assay(experiment))
  metadata <- as.data.frame(SummarizedExperiment::colData(experiment))
  X <- t(assay)
  storage.mode(X) <- "double"
  original_dimensions <- c(samples = nrow(X), taxa = ncol(X))

  sample_ids <- rownames(X)
  if (is.null(sample_ids) || is.null(rownames(metadata))) {
    stop("Qin assay and metadata must both have sample identifiers.")
  }
  metadata_index <- match(sample_ids, rownames(metadata))
  if (anyNA(metadata_index)) {
    stop("Not every Qin assay sample has a matching metadata row.")
  }
  metadata <- metadata[metadata_index, , drop = FALSE]
  if (!identical(rownames(metadata), sample_ids)) {
    stop("Qin assay/metadata alignment failed.")
  }

  colnames(X) <- make.unique(colnames(X))
  X[!is.finite(X)] <- 0
  X[X < 0] <- 0

  row_totals <- rowSums(X)
  nonempty <- is.finite(row_totals) & row_totals > 0
  X <- X[nonempty, , drop = FALSE]
  metadata <- metadata[nonempty, , drop = FALSE]
  row_totals <- row_totals[nonempty]

  converted_from_percent <- stats::median(row_totals, na.rm = TRUE) > 2
  if (converted_from_percent) {
    X <- X / 100
  }

  prevalence <- colMeans(X > 0)
  variances <- apply(X, 2L, stats::var, na.rm = TRUE)
  keep <- prevalence >= prevalence_threshold & is.finite(variances) & variances > 0
  if (!any(keep)) {
    stop("No taxa remain after prevalence and variance filtering.")
  }

  # Intentionally do not divide by rowSums(X) here. Removing taxa changes row
  # totals, and the requested analysis uses no closure adjustment.
  X <- X[, keep, drop = FALSE]

  list(
    X = X,
    metadata = metadata,
    prevalence = prevalence[keep],
    converted_from_percent = converted_from_percent,
    selected_resource = selected_resource,
    resource_candidates = resource_names %||% character(0),
    filtered_input_md5 = qin_object_md5(list(X = X, metadata = metadata)),
    original_dimensions = original_dimensions,
    filtered_dimensions = c(samples = nrow(X), taxa = ncol(X))
  )
}

make_outcome_dataset <- function(X, y, outcome, name) {
  keep_rows <- is.finite(y)
  X <- X[keep_rows, , drop = FALSE]
  y <- as.numeric(y[keep_rows])

  variances <- apply(X, 2L, stats::var, na.rm = TRUE)
  keep_columns <- is.finite(variances) & variances > 0
  X <- X[, keep_columns, drop = FALSE]

  if (nrow(X) < 20L || ncol(X) < 2L) {
    stop(name, " has insufficient observations or taxa after filtering.")
  }
  if (outcome == "binary" && !identical(sort(unique(y)), c(0, 1))) {
    stop(name, " must contain both binary classes 0 and 1.")
  }

  # No re-closure after outcome-specific variance filtering.
  list(X = X, y = y, outcome = outcome, name = name)
}

q_key <- function(q) paste0("q", sprintf("%a", as.numeric(q)))

selection_to_indices <- function(selection, feature_names) {
  if (is.null(selection) || length(selection) == 0L) {
    return(integer(0L))
  }
  if (is.logical(selection)) {
    if (length(selection) != length(feature_names)) {
      stop("A logical selection vector has the wrong length.")
    }
    return(which(selection))
  }
  if (is.character(selection)) {
    index <- match(selection, feature_names)
    return(unique(index[!is.na(index)]))
  }
  index <- unique(as.integer(selection))
  index[is.finite(index) & index >= 1L & index <= length(feature_names)]
}

lookup_q_entry <- function(x, q) {
  if (is.null(x)) return(NULL)
  keys <- c(q_key(q), as.character(q), format(q, trim = TRUE))
  if (is.list(x)) {
    for (key in keys) {
      if (!is.null(x[[key]])) return(x[[key]])
    }
  }
  NULL
}

lookup_inclusion <- function(result, method, q) {
  rates <- result$inclusion_rates %||% result$inclusion
  if (is.null(rates)) return(NULL)

  method_rates <- if (is.list(rates) && !is.null(rates[[method]])) {
    rates[[method]]
  } else {
    rates
  }
  candidate <- lookup_q_entry(method_rates, q) %||% method_rates
  if (is.numeric(candidate)) as.numeric(candidate) else NULL
}

selection_statistic <- function(result, method, q, feature_count) {
  inclusion <- lookup_inclusion(result, method, q)
  if (!is.null(inclusion) && length(inclusion) == feature_count) {
    return(inclusion)
  }

  if (method %in% c("Knockoff", "KnockoffPlus") &&
      is.numeric(result$W) && length(result$W) == feature_count) {
    return(as.numeric(result$W))
  }

  raw_p <- result$raw_pvalues %||% result$pvalues %||% result$p_values
  if (is.numeric(raw_p) && length(raw_p) == feature_count) {
    return(as.numeric(raw_p))
  }

  rep(NA_real_, feature_count)
}

result_selection_rows <- function(result, outcome_name, feature_names, q_values) {
  selections <- result$selections
  if (is.null(selections) || !is.list(selections)) {
    return(list(selected = NULL, summary = NULL, inclusion = NULL))
  }

  selected_rows <- list()
  summary_rows <- list()
  inclusion_rows <- list()
  selected_id <- summary_id <- inclusion_id <- 1L

  for (method in names(selections)) {
    for (q in q_values) {
      selected <- selection_to_indices(
        lookup_q_entry(selections[[method]], q),
        feature_names
      )
      statistic <- selection_statistic(
        result, method, q, length(feature_names)
      )

      summary_rows[[summary_id]] <- data.frame(
        Outcome = outcome_name,
        Method = method,
        q = q,
        N_Selected = length(selected),
        stringsAsFactors = FALSE
      )
      summary_id <- summary_id + 1L

      if (length(selected) > 0L) {
        selected_rows[[selected_id]] <- data.frame(
          Outcome = outcome_name,
          Method = method,
          q = q,
          Feature = feature_names[selected],
          Statistic = statistic[selected],
          stringsAsFactors = FALSE
        )
        selected_id <- selected_id + 1L
      }

      inclusion <- lookup_inclusion(result, method, q)
      if (!is.null(inclusion) && length(inclusion) == length(feature_names)) {
        inclusion_rows[[inclusion_id]] <- data.frame(
          Outcome = outcome_name,
          Method = method,
          q = q,
          Feature = feature_names,
          InclusionRate = as.numeric(inclusion),
          Selected = seq_along(feature_names) %in% selected,
          stringsAsFactors = FALSE
        )
        inclusion_id <- inclusion_id + 1L
      }
    }
  }

  list(
    selected = if (length(selected_rows)) do.call(rbind, selected_rows) else NULL,
    summary = if (length(summary_rows)) do.call(rbind, summary_rows) else NULL,
    inclusion = if (length(inclusion_rows)) do.call(rbind, inclusion_rows) else NULL
  )
}

as_diagnostic_frame <- function(x, outcome_name, source) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) {
    frame <- x
  } else if (is.list(x) && length(x) > 0L &&
             all(vapply(x, is.data.frame, logical(1L)))) {
    frame <- do.call(rbind, x)
  } else {
    return(NULL)
  }
  if (nrow(frame) == 0L) return(NULL)
  cbind(
    data.frame(Outcome = outcome_name, Source = source, stringsAsFactors = FALSE),
    frame
  )
}

safe_method_call <- function(label, expression, verbose = FALSE) {
  if (verbose) message("  Running ", label, " ...")
  tryCatch(
    list(result = force(expression), error = NULL),
    error = function(e) {
      warning(label, " failed: ", conditionMessage(e))
      list(result = NULL, error = conditionMessage(e))
    }
  )
}

run_qin_outcome <- function(dataset, config, api, seed_offset = 0L) {
  X <- dataset$X
  y <- dataset$y
  outcome <- dataset$outcome
  outcome_name <- dataset$name
  seed <- config$seed + as.integer(seed_offset)

  message(
    "Analyzing ", outcome_name, ": N=", nrow(X), ", P=", ncol(X),
    ", outcome=", outcome
  )

  negative_controls <- api$generate_gaussian_pc_controls(
    X = X,
    K = config$K,
    seed = seed + 1L,
    support = "clip",
    lower = 0,
    upper = 1,
    verbose = config$verbose
  )
  if (!is.list(negative_controls) ||
      !is.matrix(negative_controls$values) ||
      !identical(dim(negative_controls$values), dim(X))) {
    stop("Shared Gaussian-PC control generator returned an invalid object.")
  }

  xgb_control <- api$default_xgboost_control()
  if (is.list(xgb_control)) {
    xgb_control$nthread <- config$cores
  }

  method_calls <- list()
  method_calls$delporte <- safe_method_call(
    "Delporte MDS (XGBoost; Gaussian PC5, direct replacement)",
    api$run_delporte_mds(
      X = X,
      y = y,
      outcome = outcome,
      relationship = "linear",
      negative_controls = negative_controls$values,
      q_values = config$q_values,
      D = config$n_splits,
      xgb_control = xgb_control,
      K = config$K,
      support = "clip",
      seed = seed + 100L,
      verbose = config$verbose
    ),
    verbose = config$verbose
  )

  # Match the proposed method's linear predictor scale: component-wise log
  # abundance with a fixed numerical floor. This is not an ILR/CLR transform,
  # and the filtered matrix is never re-closed.
  Z_comparator <- log(pmax(X, 1e-12))
  method_calls$dai <- safe_method_call(
    "Dai MDS",
    api$run_dai_mds(
      Z = Z_comparator,
      y = y,
      outcome = outcome,
      q_values = config$q_values,
      D = config$n_splits,
      nfolds = 10L,
      seed = seed + 300L,
      verbose = config$verbose
    ),
    verbose = config$verbose
  )

  method_calls$by <- safe_method_call(
    "univariate BY",
    api$run_univariate_by(
      Z = Z_comparator,
      y = y,
      outcome = outcome,
      q_values = config$q_values
    ),
    verbose = config$verbose
  )

  method_calls$knockoffs <- safe_method_call(
    "Knockoff / Knockoff+",
    api$run_knockoffs(
      Z = Z_comparator,
      y = y,
      outcome = outcome,
      q_values = config$q_values,
      seed = seed + 400L,
      nfolds = 10L
    ),
    verbose = config$verbose
  )

  # Some shared comparators return a valid empty result together with an
  # internal failure message (notably second-order knockoffs when covariance
  # estimation is ill-conditioned). Preserve that message in method_errors.csv.
  for (method_name in names(method_calls)) {
    internal_error <- method_calls[[method_name]]$result$error %||% NA_character_
    if (length(internal_error) == 1L && !is.na(internal_error) &&
        nzchar(internal_error)) {
      method_calls[[method_name]]$error <- internal_error
    }
  }

  list(
    dataset = dataset,
    negative_controls = negative_controls,
    methods = lapply(method_calls, `[[`, "result"),
    errors = vapply(
      method_calls,
      function(x) x$error %||% NA_character_,
      character(1L)
    )
  )
}

combine_non_null_frames <- function(frames) {
  frames <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, frames)
  if (length(frames) == 0L) return(data.frame())
  all_names <- unique(unlist(lapply(frames, names), use.names = FALSE))
  frames <- lapply(frames, function(x) {
    missing <- setdiff(all_names, names(x))
    for (name in missing) x[[name]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, frames)
}

extract_by_table <- function(result, outcome_name, feature_names) {
  if (is.null(result)) return(NULL)
  raw <- result$raw_pvalues %||% result$pvalues %||% result$p_values %||%
    result$p_value
  adjusted <- result$adjusted_pvalues %||% result$by_pvalues %||%
    result$adjusted_p_values %||% result$pvalues_adjusted

  if (!is.numeric(raw) || length(raw) != length(feature_names)) return(NULL)
  if (!is.numeric(adjusted) || length(adjusted) != length(feature_names)) {
    adjusted <- stats::p.adjust(raw, method = "BY")
  }

  data.frame(
    Outcome = outcome_name,
    Feature = feature_names,
    PValue = as.numeric(raw),
    BY_AdjustedPValue = as.numeric(adjusted),
    stringsAsFactors = FALSE
  )
}

write_qin_outputs <- function(results, config, input_info) {
  dir.create(config$output, recursive = TRUE, showWarnings = FALSE)

  selected_frames <- list()
  summary_frames <- list()
  inclusion_frames <- list()
  diagnostic_frames <- list()
  nc_diagnostic_frames <- list()
  by_frames <- list()
  errors <- list()
  status_frames <- list()
  id <- 1L

  for (outcome_name in names(results)) {
    outcome_result <- results[[outcome_name]]
    feature_names <- colnames(outcome_result$dataset$X)

    nc_diag <- outcome_result$negative_controls$diagnostics
    if (is.data.frame(nc_diag) && nrow(nc_diag) > 0L) {
      nc_diagnostic_frames[[length(nc_diagnostic_frames) + 1L]] <- cbind(
        data.frame(Outcome = outcome_name, stringsAsFactors = FALSE),
        nc_diag
      )
    }

    for (source in names(outcome_result$methods)) {
      method_result <- outcome_result$methods[[source]]
      if (is.null(method_result)) next

      diagnostic_frames[[id]] <- as_diagnostic_frame(
        method_result$split_diagnostics %||% method_result$mirror_diagnostics,
        outcome_name,
        source
      )
      error_text <- outcome_result$errors[[source]] %||% NA_character_
      method_failed <- length(error_text) == 1L && !is.na(error_text) &&
        nzchar(trimws(error_text))
      if (method_failed) {
        id <- id + 1L
        next
      }

      tables <- result_selection_rows(
        method_result,
        outcome_name,
        feature_names,
        config$q_values
      )
      selected_frames[[id]] <- tables$selected
      summary_frames[[id]] <- tables$summary
      inclusion_frames[[id]] <- tables$inclusion
      if (identical(source, "by")) {
        by_frames[[length(by_frames) + 1L]] <- extract_by_table(
          method_result, outcome_name, feature_names
        )
      }
      id <- id + 1L
    }

    failed <- names(outcome_result$errors)[!is.na(outcome_result$errors)]
    method_names <- names(outcome_result$errors)
    if (length(method_names) > 0L) {
      error_values <- unname(outcome_result$errors)
      status_frames[[outcome_name]] <- data.frame(
        Outcome = outcome_name,
        MethodCall = method_names,
        Status = ifelse(is.na(error_values), "success", "failed"),
        Error = error_values,
        stringsAsFactors = FALSE
      )
    }
    if (length(failed) > 0L) {
      errors[[outcome_name]] <- data.frame(
        Outcome = outcome_name,
        MethodCall = failed,
        Error = unname(outcome_result$errors[failed]),
        stringsAsFactors = FALSE
      )
    }
  }

  selected <- combine_non_null_frames(selected_frames)
  summary <- combine_non_null_frames(summary_frames)
  inclusion <- combine_non_null_frames(inclusion_frames)
  diagnostics <- combine_non_null_frames(diagnostic_frames)
  nc_diagnostics <- combine_non_null_frames(nc_diagnostic_frames)
  by_pvalues <- combine_non_null_frames(by_frames)
  error_table <- combine_non_null_frames(errors)
  status_table <- combine_non_null_frames(status_frames)

  if (ncol(selected) == 0L) {
    selected <- data.frame(
      Outcome = character(0), Method = character(0), q = numeric(0),
      Feature = character(0), Statistic = numeric(0)
    )
  }
  if (ncol(summary) == 0L) {
    summary <- data.frame(
      Outcome = character(0), Method = character(0), q = numeric(0),
      N_Selected = integer(0)
    )
  }
  if (ncol(inclusion) == 0L) {
    inclusion <- data.frame(
      Outcome = character(0), Method = character(0), q = numeric(0),
      Feature = character(0), InclusionRate = numeric(0),
      Selected = logical(0)
    )
  }
  if (ncol(diagnostics) == 0L) {
    diagnostics <- data.frame(Outcome = character(0), Source = character(0))
  }
  if (ncol(nc_diagnostics) == 0L) {
    nc_diagnostics <- data.frame(
      Outcome = character(0), feature = integer(0), heldout_fold = integer(0),
      n_generated = integer(0), K_effective = integer(0),
      residual_sd = numeric(0), n_below_support = integer(0),
      n_above_support = integer(0), n_clipped = integer(0),
      clip_rate = numeric(0), n_at_or_below_log_floor = integer(0),
      used_fallback = logical(0)
    )
  }
  if (ncol(by_pvalues) == 0L) {
    by_pvalues <- data.frame(
      Outcome = character(0), Feature = character(0), PValue = numeric(0),
      BY_AdjustedPValue = numeric(0)
    )
  }
  if (ncol(error_table) == 0L) {
    error_table <- data.frame(
      Outcome = character(0), MethodCall = character(0), Error = character(0)
    )
  }
  if (ncol(status_table) == 0L) {
    status_table <- data.frame(
      Outcome = character(0), MethodCall = character(0),
      Status = character(0), Error = character(0)
    )
  }

  utils::write.csv(
    selected,
    file.path(config$output, "selected_taxa.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summary,
    file.path(config$output, "selection_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    inclusion,
    file.path(config$output, "inclusion_rates.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    diagnostics,
    file.path(config$output, "split_diagnostics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    nc_diagnostics,
    file.path(config$output, "negative_control_diagnostics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    by_pvalues,
    file.path(config$output, "BY_pvalues.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    error_table,
    file.path(config$output, "method_errors.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    status_table,
    file.path(config$output, "method_status.csv"),
    row.names = FALSE
  )

  session <- utils::sessionInfo()
  writeLines(
    capture.output(session),
    file.path(config$output, "session_info.txt")
  )

  package_names <- c(
    "curatedMetagenomicData", "SummarizedExperiment", "glmnet",
    "xgboost", "knockoff"
  )
  package_versions <- vapply(
    package_names,
    function(package) as.character(utils::packageVersion(package)),
    character(1L)
  )
  manifest <- data.frame(
    Key = c(
      "run_timestamp_utc", "dataset", "curated_resource_selected",
      "curated_resource_candidates", "filtered_input_md5", "PMID", "counts",
      "analysis_features",
      "input_scale", "prevalence_threshold", "post_filter_closure",
      "negative_control_model", "negative_control_pcs",
      "negative_control_transform", "perturbation", "q_values", "n_splits",
      "seed", "cores", "shared_workflow_version", "shared_workflow_code_hash",
      "qin_workflow_version", "qin_workflow_code_hash", "original_dimensions",
      "filtered_dimensions", paste0("package_", package_names)
    ),
    Value = c(
      format(Sys.time(), tz = "UTC", usetz = TRUE),
      "QinN_2014.relative_abundance",
      input_info$selected_resource,
      paste(input_info$resource_candidates, collapse = " | "),
      input_info$filtered_input_md5,
      "25079328",
      "FALSE",
      "taxa_only",
      "relative_abundance_proportion",
      format(config$prevalence, scientific = FALSE),
      "none",
      "cross_fitted_Gaussian_regression",
      as.character(config$K),
      "first_K_PCs_of_raw_X_minus_j_no_ILR_or_CLR",
      "replace_j_only_no_closure_adjustment",
      paste(config$q_values, collapse = ","),
      as.character(config$n_splits),
      as.character(config$seed),
      as.character(config$cores),
      input_info$shared_workflow_version,
      as.character(input_info$shared_workflow_code_hash),
      QIN_WORKFLOW_VERSION,
      unname(tools::md5sum(QIN_SCRIPT_FILE)),
      paste(input_info$original_dimensions, collapse = " x "),
      paste(input_info$filtered_dimensions, collapse = " x "),
      unname(package_versions)
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    manifest,
    file.path(config$output, "run_manifest.csv"),
    row.names = FALSE
  )

  saveRDS(
    list(
      config = config,
      input = input_info,
      outcomes = results,
      selected_taxa = selected,
      selection_summary = summary,
      inclusion_rates = inclusion,
      split_diagnostics = diagnostics,
      negative_control_diagnostics = nc_diagnostics,
      BY_pvalues = by_pvalues,
      method_errors = error_table,
      method_status = status_table,
      manifest = manifest,
      session_info = session
    ),
    file.path(config$output, "qin_2014_full_results.rds")
  )

  invisible(list(
    selected_taxa = selected,
    selection_summary = summary,
    inclusion_rates = inclusion,
    errors = error_table,
    method_status = status_table,
    manifest = manifest
  ))
}

run_qin_2014_analysis <- function(config = parse_qin_cli(character())) {
  check_qin_packages()
  api <- load_shared_simulation_api()

  if (isTRUE(config$dry_run)) {
    message("Dry run successful; no data were downloaded and no results were written.")
    message("Repository root: ", QIN_REPO_ROOT)
    message("Shared methods: ", QIN_SHARED_SCRIPT)
    message("Output directory: ", config$output)
    message("q: ", paste(config$q_values, collapse = ", "))
    message("splits: ", config$n_splits, "; cores: ", config$cores)
    message("negative controls: Gaussian raw-PC K=", config$K,
            ", direct replacement, no closure")
    return(invisible(config))
  }

  input <- load_qin_2014(config$prevalence)
  input$shared_workflow_version <- api$WORKFLOW_VERSION
  input$shared_workflow_code_hash <- api$workflow_code_hash()
  albumin <- make_outcome_dataset(
    input$X,
    extract_albumin(input$metadata),
    outcome = "continuous",
    name = "Albumin"
  )
  cirrhosis <- make_outcome_dataset(
    input$X,
    extract_cirrhosis(input$metadata),
    outcome = "binary",
    name = "Cirrhosis"
  )

  message(
    "Qin data prepared without post-filter closure: N=", nrow(input$X),
    ", P=", ncol(input$X)
  )
  results <- list(
    Albumin = run_qin_outcome(albumin, config, api, seed_offset = 10000L),
    Cirrhosis = run_qin_outcome(cirrhosis, config, api, seed_offset = 20000L)
  )

  output_tables <- write_qin_outputs(results, config, input)
  message("QinN_2014 analysis complete. Results written to: ", config$output)
  invisible(list(results = results, tables = output_tables))
}

qin_main <- function() {
  config <- parse_qin_cli()
  if (isTRUE(config$help)) {
    cat(qin_usage(), "\n")
    return(invisible(NULL))
  }
  run_qin_2014_analysis(config)
}

if (qin_is_direct_execution()) {
  qin_main()
}
