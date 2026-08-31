#!/usr/bin/env Rscript

# Summarize and plot the unified simulation output.
#
# Command-line examples (paths are relative to the repository root):
#   Rscript reproducibility/simulation/plot_simulation_results.R
#   Rscript reproducibility/simulation/plot_simulation_results.R \
#     --input reproducibility/results/main/single_dirichlet/replicate_metrics.csv \
#     --output reproducibility/results/main/single_dirichlet/figures \
#     --study main --q 0.05,0.10
#
# Sourcing this file defines functions only; it does not read data or draw plots.

REQUIRED_METRIC_COLUMNS <- c(
  "study",
  "generator",
  "J",
  "n",
  "n_signal",
  "outcome",
  "relationship",
  "rho",
  "replicate",
  "method",
  "q",
  "fdp",
  "power",
  "n_selected"
)

.plot_script_file <- local({
  frames <- rev(sys.frames())
  for (frame in frames) {
    if (!is.null(frame$ofile)) {
      return(normalizePath(frame$ofile, mustWork = FALSE))
    }
  }

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

  file_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  if (length(file_argument) > 0L) {
    return(normalizePath(
      sub("^--file=", "", file_argument[[1L]]),
      mustWork = FALSE
    ))
  }

  NA_character_
})

project_paths <- function(script_file = .plot_script_file) {
  if (length(script_file) != 1L || is.na(script_file)) {
    stop(
      "Could not determine the plotting script location. ",
      "Run this file with Rscript or supply absolute input/output paths.",
      call. = FALSE
    )
  }

  script_directory <- dirname(script_file)
  reproducibility_root <- dirname(script_directory)
  repository_root <- dirname(reproducibility_root)

  list(
    repository_root = normalizePath(repository_root, mustWork = TRUE),
    reproducibility_root = normalizePath(
      reproducibility_root,
      mustWork = TRUE
    )
  )
}

is_absolute_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
}

resolve_repository_path <- function(path, repository_root) {
  path <- path.expand(path)
  if (!is_absolute_path(path)) {
    path <- file.path(repository_root, path)
  }
  normalizePath(path, mustWork = FALSE)
}

split_option_values <- function(values) {
  if (length(values) == 0L) return(character(0))
  pieces <- unlist(strsplit(values, ",", fixed = TRUE), use.names = FALSE)
  pieces <- trimws(pieces)
  unique(pieces[nzchar(pieces)])
}

parse_command_line <- function(args) {
  options <- list(
    input = NULL,
    output = NULL,
    study = character(0),
    q = character(0),
    allow_failures = FALSE,
    help = FALSE
  )

  take_value <- function(index, option_name) {
    if (index >= length(args)) {
      stop("Missing value after ", option_name, ".", call. = FALSE)
    }
    args[[index + 1L]]
  }

  i <- 1L
  while (i <= length(args)) {
    argument <- args[[i]]

    if (argument %in% c("--help", "-h")) {
      options$help <- TRUE
    } else if (identical(argument, "--allow-failures")) {
      options$allow_failures <- TRUE
    } else if (startsWith(argument, "--input=")) {
      options$input <- sub("^--input=", "", argument)
    } else if (identical(argument, "--input")) {
      options$input <- take_value(i, "--input")
      i <- i + 1L
    } else if (startsWith(argument, "--output=")) {
      options$output <- sub("^--output=", "", argument)
    } else if (identical(argument, "--output")) {
      options$output <- take_value(i, "--output")
      i <- i + 1L
    } else if (startsWith(argument, "--study=")) {
      options$study <- c(
        options$study,
        sub("^--study=", "", argument)
      )
    } else if (identical(argument, "--study")) {
      options$study <- c(options$study, take_value(i, "--study"))
      i <- i + 1L
    } else if (startsWith(argument, "--q=")) {
      options$q <- c(options$q, sub("^--q=", "", argument))
    } else if (identical(argument, "--q")) {
      options$q <- c(options$q, take_value(i, "--q"))
      i <- i + 1L
    } else {
      stop(
        "Unknown argument: ", argument,
        ". Use --help to see the supported interface.",
        call. = FALSE
      )
    }

    i <- i + 1L
  }

  options$study <- split_option_values(options$study)
  options$q <- split_option_values(options$q)
  options
}

print_usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript reproducibility/simulation/plot_simulation_results.R [options]\n\n",
      "Options:\n",
      "  --input PATH    Unified replicate-level CSV (required).\n",
      "  --output DIR    Figure/output directory. Default: figures beside input.\n",
      "  --study VALUE   Keep one or more study values (repeat the option or\n",
      "                  provide a comma-separated list).\n",
      "  --q VALUE       Keep one or more FDR targets (repeat the option or\n",
      "                  provide a comma-separated list).\n",
      "  --allow-failures  Plot finite replicates even if method failures exist.\n",
      "  --help, -h      Show this help.\n\n",
      "Relative paths are resolved from the repository root.\n"
    )
  )
}

coerce_numeric_column <- function(x, column_name) {
  original <- trimws(as.character(x))
  is_missing <- is.na(x) | !nzchar(original) |
    original %in% c("NA", "NaN", "na", "nan")

  converted <- suppressWarnings(as.numeric(original))
  invalid <- !is_missing & is.na(converted)
  if (any(invalid)) {
    bad_rows <- head(which(invalid), 5L)
    stop(
      "Column '", column_name, "' contains non-numeric values at row(s): ",
      paste(bad_rows, collapse = ", "),
      call. = FALSE
    )
  }

  converted[is_missing] <- NA_real_
  converted
}

validate_replicate_metrics <- function(metrics) {
  missing_columns <- setdiff(REQUIRED_METRIC_COLUMNS, names(metrics))
  if (length(missing_columns) > 0L) {
    stop(
      "Input is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  metrics <- metrics[, REQUIRED_METRIC_COLUMNS, drop = FALSE]
  if (nrow(metrics) == 0L) {
    stop("The input CSV contains no data rows.", call. = FALSE)
  }

  text_columns <- c(
    "study", "generator", "outcome", "relationship", "method"
  )
  for (column_name in text_columns) {
    values <- trimws(as.character(metrics[[column_name]]))
    values[!nzchar(values)] <- NA_character_
    metrics[[column_name]] <- values
  }

  for (column_name in c("study", "method")) {
    if (anyNA(metrics[[column_name]])) {
      bad_rows <- head(which(is.na(metrics[[column_name]])), 5L)
      stop(
        "Column '", column_name, "' is blank at row(s): ",
        paste(bad_rows, collapse = ", "),
        call. = FALSE
      )
    }
  }

  numeric_columns <- c(
    "J", "n", "n_signal", "rho", "replicate", "q",
    "fdp", "power", "n_selected"
  )
  for (column_name in numeric_columns) {
    metrics[[column_name]] <- coerce_numeric_column(
      metrics[[column_name]],
      column_name
    )
  }

  assert_finite <- function(column_name, allow_na = FALSE) {
    values <- metrics[[column_name]]
    invalid <- !is.finite(values)
    if (allow_na) invalid <- !is.na(values) & invalid
    if (any(invalid)) {
      stop(
        "Column '", column_name, "' contains missing or infinite values.",
        call. = FALSE
      )
    }
  }

  for (column_name in c("J", "n", "replicate", "q")) {
    assert_finite(column_name)
  }
  for (column_name in c("n_signal", "rho", "fdp", "power", "n_selected")) {
    assert_finite(column_name, allow_na = TRUE)
  }

  if (any(metrics$J <= 0 | metrics$n <= 0 | metrics$replicate <= 0)) {
    stop("J, n, and replicate must be positive.", call. = FALSE)
  }
  if (any(abs(metrics$J - round(metrics$J)) > 1e-8) ||
      any(abs(metrics$n - round(metrics$n)) > 1e-8) ||
      any(abs(metrics$replicate - round(metrics$replicate)) > 1e-8)) {
    stop("J, n, and replicate must be integer-valued.", call. = FALSE)
  }
  if (any(!is.na(metrics$n_signal) & metrics$n_signal < 0) ||
      any(!is.na(metrics$n_selected) & metrics$n_selected < 0)) {
    stop("n_signal and n_selected cannot be negative.", call. = FALSE)
  }
  if (any(metrics$q <= 0 | metrics$q >= 1)) {
    stop("All q values must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (any(!is.na(metrics$fdp) &
          (metrics$fdp < -1e-10 | metrics$fdp > 1 + 1e-10)) ||
      any(!is.na(metrics$power) &
          (metrics$power < -1e-10 | metrics$power > 1 + 1e-10))) {
    stop("FDP and power values must be between 0 and 1.", call. = FALSE)
  }
  if (any(!is.na(metrics$rho) & abs(metrics$rho) > 1 + 1e-10)) {
    stop("Finite rho values must be between -1 and 1.", call. = FALSE)
  }

  identity_columns <- c(
    "study", "generator", "J", "n", "n_signal", "outcome",
    "relationship", "rho", "replicate", "method", "q"
  )
  identity_parts <- lapply(metrics[identity_columns], function(x) {
    ifelse(is.na(x), "<NA>", as.character(x))
  })
  identity_key <- do.call(
    paste,
    c(identity_parts, list(sep = "\034"))
  )
  duplicate_rows <- which(duplicated(identity_key))
  if (length(duplicate_rows) > 0L) {
    stop(
      "Duplicate scenario/replicate/method/q rows were found; first duplicate ",
      "is at input row ", duplicate_rows[[1L]], ".",
      call. = FALSE
    )
  }

  metrics$J <- as.integer(round(metrics$J))
  metrics$n <- as.integer(round(metrics$n))
  metrics$replicate <- as.integer(round(metrics$replicate))
  metrics$fdp <- pmin(pmax(metrics$fdp, 0), 1)
  metrics$power <- pmin(pmax(metrics$power, 0), 1)
  metrics
}

read_replicate_metrics <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }

  metrics <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  )
  validate_replicate_metrics(metrics)
}

filter_replicate_metrics <- function(metrics, studies = NULL, q_values = NULL) {
  keep <- rep(TRUE, nrow(metrics))

  if (!is.null(studies) && length(studies) > 0L) {
    studies <- split_option_values(studies)
    available <- unique(metrics$study)
    requested_key <- tolower(studies)
    available_key <- tolower(available)
    unmatched <- studies[!requested_key %in% available_key]
    if (length(unmatched) > 0L) {
      stop(
        "Requested --study value(s) not found: ",
        paste(unmatched, collapse = ", "),
        ". Available values: ", paste(sort(available), collapse = ", "),
        call. = FALSE
      )
    }
    keep <- keep & tolower(metrics$study) %in% requested_key
  }

  if (!is.null(q_values) && length(q_values) > 0L) {
    q_text <- split_option_values(as.character(q_values))
    requested_q <- suppressWarnings(as.numeric(q_text))
    if (anyNA(requested_q) || any(requested_q <= 0 | requested_q >= 1)) {
      stop("Every --q value must be numeric and between 0 and 1.", call. = FALSE)
    }

    q_match <- vapply(requested_q, function(target) {
      any(abs(metrics$q - target) <= 1e-10)
    }, logical(1))
    if (any(!q_match)) {
      stop(
        "Requested --q value(s) not found: ",
        paste(q_text[!q_match], collapse = ", "),
        ". Available values: ",
        paste(sort(unique(metrics$q)), collapse = ", "),
        call. = FALSE
      )
    }

    keep_q <- Reduce(
      `|`,
      lapply(requested_q, function(target) {
        abs(metrics$q - target) <= 1e-10
      })
    )
    keep <- keep & keep_q
  }

  filtered <- metrics[keep, , drop = FALSE]
  if (nrow(filtered) == 0L) {
    stop("No input rows remain after applying filters.", call. = FALSE)
  }
  rownames(filtered) <- NULL
  filtered
}

summarize_vector <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) {
    return(c(
      n = 0,
      mean = NA_real_,
      median = NA_real_,
      mcse = NA_real_,
      q025 = NA_real_,
      q975 = NA_real_
    ))
  }

  quantiles <- stats::quantile(
    x,
    probs = c(0.025, 0.975),
    names = FALSE,
    type = 7
  )
  c(
    n = n,
    mean = mean(x),
    median = stats::median(x),
    mcse = if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_,
    q025 = quantiles[[1L]],
    q975 = quantiles[[2L]]
  )
}

summarize_replicate_metrics <- function(metrics) {
  group_columns <- c(
    "study", "generator", "J", "n", "n_signal", "outcome",
    "relationship", "rho", "method", "q"
  )
  key_parts <- lapply(metrics[group_columns], function(x) {
    ifelse(is.na(x), "<NA>", as.character(x))
  })
  group_key <- do.call(paste, c(key_parts, list(sep = "\034")))
  groups <- split(seq_len(nrow(metrics)), group_key, drop = TRUE)

  rows <- lapply(groups, function(index) {
    fdp_summary <- summarize_vector(metrics$fdp[index])
    power_summary <- summarize_vector(metrics$power[index])
    group_values <- metrics[index[[1L]], group_columns, drop = FALSE]

    cbind(
      group_values,
      data.frame(
        n_rows = length(index),
        n_replicates = length(unique(metrics$replicate[index])),
        n_fdp = as.integer(fdp_summary[["n"]]),
        mean_fdp = fdp_summary[["mean"]],
        median_fdp = fdp_summary[["median"]],
        mcse_fdp = fdp_summary[["mcse"]],
        q025_fdp = fdp_summary[["q025"]],
        q975_fdp = fdp_summary[["q975"]],
        n_power = as.integer(power_summary[["n"]]),
        mean_power = power_summary[["mean"]],
        median_power = power_summary[["median"]],
        mcse_power = power_summary[["mcse"]],
        q025_power = power_summary[["q025"]],
        q975_power = power_summary[["q975"]],
        check.names = FALSE
      )
    )
  })

  summary_table <- do.call(rbind, rows)
  rownames(summary_table) <- NULL
  order_columns <- c(
    "study", "generator", "J", "n", "n_signal", "outcome",
    "relationship", "rho", "q", "method"
  )
  # Drop the column names before do.call(): otherwise the grouping column
  # named "method" is matched to order()'s own method argument.
  order_arguments <- unname(summary_table[order_columns])
  ordering <- do.call(
    order,
    c(order_arguments, list(na.last = TRUE))
  )
  summary_table[ordering, , drop = FALSE]
}

pretty_label <- function(x, missing = "unspecified") {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- missing
  gsub("_", " ", x, fixed = TRUE)
}

format_number <- function(x, digits = 4L) {
  out <- format(x, trim = TRUE, scientific = FALSE, digits = digits)
  out[is.na(x)] <- "unspecified"
  out
}

make_main_panel_labels <- function(data) {
  paste0(
    pretty_label(data$study), " | ", pretty_label(data$generator), "\n",
    "J=", data$J,
    ", n=", data$n,
    ", ", pretty_label(data$outcome),
    ", ", pretty_label(data$relationship),
    ", q=", format_number(data$q)
  )
}

make_correlation_panel_labels <- function(data) {
  paste0(
    pretty_label(data$study), " | ", pretty_label(data$generator), "\n",
    "J=", data$J,
    ", n=", data$n,
    ", signals=", format_number(data$n_signal),
    ", ", pretty_label(data$outcome),
    ", ", pretty_label(data$relationship),
    ", q=", format_number(data$q)
  )
}

classify_plot_rows <- function(data) {
  labels <- tolower(paste(
    ifelse(is.na(data$study), "", data$study),
    ifelse(is.na(data$generator), "", data$generator)
  ))
  correlation_hint <- grepl(
    "correlation|latent[ _-]*block|block[ _-]*gaussian|logistic[ _-]*normal",
    labels
  )
  is_correlation <- is.finite(data$rho) | correlation_hint
  outcome_key <- tolower(trimws(data$outcome))

  list(
    main = !is_correlation &
      data$J %in% c(100L, 1000L) &
      outcome_key %in% c("continuous", "binary"),
    correlation = is_correlation & data$J == 100L & is.finite(data$rho)
  )
}

metric_plot_data <- function(data, metric, plot_type) {
  stopifnot(metric %in% c("fdp", "power"))
  stopifnot(plot_type %in% c("main", "correlation"))

  mean_column <- paste0("mean_", metric)
  mcse_column <- paste0("mcse_", metric)
  data$estimate <- data[[mean_column]]
  data$mcse <- data[[mcse_column]]
  data$ci_low <- pmax(0, data$estimate - 1.96 * data$mcse)
  data$ci_high <- pmin(1, data$estimate + 1.96 * data$mcse)
  one_observation <- is.finite(data$estimate) & !is.finite(data$mcse)
  data$ci_low[one_observation] <- data$estimate[one_observation]
  data$ci_high[one_observation] <- data$estimate[one_observation]

  if (plot_type == "main") {
    data$panel <- make_main_panel_labels(data)
    signal_label <- ifelse(
      is.na(data$n_signal),
      "unspecified",
      format_number(data$n_signal)
    )
    numeric_signals <- sort(unique(data$n_signal[is.finite(data$n_signal)]))
    signal_levels <- c(format_number(numeric_signals))
    if (anyNA(data$n_signal)) signal_levels <- c(signal_levels, "unspecified")
    data$x_value <- factor(signal_label, levels = unique(signal_levels))
    x_title <- "Number of signal features"
  } else {
    data$panel <- make_correlation_panel_labels(data)
    rho_label <- ifelse(
      is.na(data$rho),
      "unspecified",
      format_number(data$rho)
    )
    numeric_rho <- sort(unique(data$rho[is.finite(data$rho)]))
    rho_levels <- c(format_number(numeric_rho))
    if (anyNA(data$rho)) rho_levels <- c(rho_levels, "unspecified")
    data$x_value <- factor(rho_label, levels = unique(rho_levels))
    x_title <- "Latent within-block correlation (rho)"
  }

  method_levels <- sort(unique(as.character(data$method)))
  data$method <- factor(data$method, levels = method_levels)

  list(data = data, x_title = x_title, method_levels = method_levels)
}

build_metric_plot <- function(data, metric, plot_type) {
  prepared <- metric_plot_data(data, metric, plot_type)
  plot_data <- prepared$data
  finite_data <- plot_data[is.finite(plot_data$estimate), , drop = FALSE]
  if (nrow(finite_data) == 0L) {
    warning("No finite ", metric, " summaries are available for ", plot_type, ".")
    return(NULL)
  }

  method_levels <- prepared$method_levels
  shape_values <- setNames(
    rep(c(16, 17, 15, 18, 8, 3, 4, 7, 9, 10, 11, 12),
        length.out = length(method_levels)),
    method_levels
  )
  linetype_values <- setNames(
    rep(c("solid", "dashed", "dotdash", "longdash", "twodash", "dotted"),
        length.out = length(method_levels)),
    method_levels
  )
  dodge <- ggplot2::position_dodge(width = 0.45)

  plot <- ggplot2::ggplot(
    finite_data,
    ggplot2::aes(
      x = x_value,
      y = estimate,
      colour = method,
      shape = method,
      linetype = method,
      group = method
    )
  )

  line_key <- paste(finite_data$panel, finite_data$method, sep = "\034")
  x_per_line <- tapply(
    as.character(finite_data$x_value),
    line_key,
    function(x) length(unique(x))
  )
  line_data <- finite_data[line_key %in% names(x_per_line[x_per_line > 1L]), ,
                           drop = FALSE]
  if (nrow(line_data) > 0L) {
    plot <- plot + ggplot2::geom_line(
      data = line_data,
      position = dodge,
      linewidth = 0.45,
      alpha = 0.85,
      na.rm = TRUE
    )
  }

  plot <- plot +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      position = dodge,
      width = 0.18,
      linewidth = 0.45,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      position = dodge,
      size = 2.4,
      stroke = 0.8,
      na.rm = TRUE
    )

  if (metric == "fdp") {
    target_lines <- unique(finite_data[c("panel", "q")])
    plot <- plot + ggplot2::geom_hline(
      data = target_lines,
      ggplot2::aes(yintercept = q),
      inherit.aes = FALSE,
      colour = "#B2182B",
      linetype = "dashed",
      linewidth = 0.55
    )
  }

  y_title <- if (metric == "fdp") "False discovery proportion" else "Power"
  title_prefix <- if (plot_type == "main") {
    "Dirichlet simulations"
  } else {
    "J = 100 correlation simulations"
  }

  plot +
    ggplot2::facet_wrap(~panel, scales = "free_x") +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = ggplot2::expansion(mult = c(0.02, 0.05))
    ) +
    ggplot2::scale_shape_manual(
      values = shape_values,
      labels = function(x) pretty_label(x)
    ) +
    ggplot2::scale_linetype_manual(
      values = linetype_values,
      labels = function(x) pretty_label(x)
    ) +
    ggplot2::scale_colour_discrete(labels = function(x) pretty_label(x)) +
    ggplot2::labs(
      title = paste(title_prefix, y_title, sep = ": "),
      x = prepared$x_title,
      y = y_title,
      colour = NULL,
      shape = NULL,
      linetype = NULL,
      caption = paste0(
        "Points are Monte Carlo means; error bars are mean +/- 1.96 MCSE. ",
        "Medians and empirical 2.5th/97.5th percentiles are in plot_summary.csv.",
        " Any explicitly allowed failed rows are excluded; counts are recorded there.",
        if (metric == "fdp") " Dashed red lines mark the target q." else ""
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey75"),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8),
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    )
}

save_plot_pair <- function(plot, output_directory, stem, n_panels) {
  if (is.null(plot)) return(character(0))

  panel_columns <- if (n_panels <= 1L) 1L else if (n_panels <= 4L) 2L else 3L
  panel_rows <- ceiling(n_panels / panel_columns)
  width <- max(7.5, 5.2 * panel_columns)
  height <- max(5.2, 3.9 * panel_rows + 1.2)

  png_path <- file.path(output_directory, paste0(stem, ".png"))
  pdf_path <- file.path(output_directory, paste0(stem, ".pdf"))
  ggplot2::ggsave(
    png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  ggplot2::ggsave(
    pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
  c(png_path, pdf_path)
}

plot_simulation_results <- function(
    input,
    output,
    studies = NULL,
    q_values = NULL,
    allow_failures = FALSE,
    repository_root = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "The ggplot2 package is required. Install it before running this script.",
      call. = FALSE
    )
  }

  if (is.null(repository_root)) {
    repository_root <- project_paths()$repository_root
  }
  input_path <- resolve_repository_path(input, repository_root)
  output_directory <- resolve_repository_path(output, repository_root)

  metrics <- read_replicate_metrics(input_path)
  metrics <- filter_replicate_metrics(
    metrics,
    studies = studies,
    q_values = q_values
  )
  failed_rows <- !is.finite(metrics$fdp) | !is.finite(metrics$power) |
    !is.finite(metrics$n_selected)
  if (any(failed_rows) && !isTRUE(allow_failures)) {
    stop(
      sum(failed_rows),
      " retained row(s) have failed-method metrics. Inspect ",
      "method_diagnostics.csv and resolve the failures, or rerun with ",
      "--allow-failures to plot finite replicates explicitly.",
      call. = FALSE
    )
  }
  if (any(failed_rows)) {
    warning(
      "Plotting finite replicates after excluding ", sum(failed_rows),
      " failed-method row(s). Counts are recorded in plot_summary.csv.",
      call. = FALSE
    )
  }
  summary_table <- summarize_replicate_metrics(metrics)

  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_directory)) {
    stop("Could not create output directory: ", output_directory, call. = FALSE)
  }

  summary_path <- file.path(output_directory, "plot_summary.csv")
  utils::write.csv(
    summary_table,
    summary_path,
    row.names = FALSE,
    na = ""
  )

  row_classes <- classify_plot_rows(summary_table)
  main_data <- summary_table[row_classes$main, , drop = FALSE]
  correlation_data <- summary_table[row_classes$correlation, , drop = FALSE]
  generated_files <- summary_path

  if (nrow(main_data) > 0L) {
    n_panels <- length(unique(make_main_panel_labels(main_data)))
    for (metric in c("fdp", "power")) {
      plot <- build_metric_plot(main_data, metric, "main")
      generated_files <- c(
        generated_files,
        save_plot_pair(
          plot,
          output_directory,
          paste0("main_", metric),
          n_panels
        )
      )
    }
  }

  if (nrow(correlation_data) > 0L) {
    n_panels <- length(unique(make_correlation_panel_labels(correlation_data)))
    for (metric in c("fdp", "power")) {
      plot <- build_metric_plot(correlation_data, metric, "correlation")
      generated_files <- c(
        generated_files,
        save_plot_pair(
          plot,
          output_directory,
          paste0("correlation_", metric),
          n_panels
        )
      )
    }
  }

  if (nrow(main_data) == 0L && nrow(correlation_data) == 0L) {
    stop(
      "No plottable rows were found. Main plots require J in {100, 1000}, ",
      "continuous/binary outcomes, and missing rho; correlation plots require ",
      "J = 100 and a finite rho.",
      call. = FALSE
    )
  }

  ignored_rows <- nrow(summary_table) - nrow(main_data) - nrow(correlation_data)
  if (ignored_rows > 0L) {
    warning(
      ignored_rows,
      " summary row(s) did not match the requested main or correlation panels."
    )
  }

  message("Read: ", input_path)
  message(
    "Retained ", nrow(metrics), " replicate-level row(s) across ",
    length(unique(metrics$study)), " study value(s)."
  )
  message("Wrote: ", summary_path)
  for (path in setdiff(generated_files, summary_path)) {
    message("Wrote: ", path)
  }

  invisible(list(
    input = input_path,
    output_directory = output_directory,
    summary = summary_table,
    files = generated_files
  ))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- parse_command_line(args)
  if (isTRUE(options$help)) {
    print_usage()
    return(invisible(NULL))
  }

  paths <- project_paths()
  if (is.null(options$input)) {
    stop("--input is required. Use --help for an example.", call. = FALSE)
  }
  input <- options$input
  output <- if (is.null(options$output)) {
    file.path(dirname(input), "figures")
  } else {
    options$output
  }

  plot_simulation_results(
    input = input,
    output = output,
    studies = options$study,
    q_values = options$q,
    allow_failures = options$allow_failures,
    repository_root = paths$repository_root
  )
}

if (sys.nframe() == 0L) {
  main()
}
