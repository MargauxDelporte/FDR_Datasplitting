# Reproducible compositional FDR simulations
#
# This file is intentionally self-contained.  It contains the proposed method,
# all comparison methods, both abundance generators, checkpointing, and the
# command-line runner used for the simulation studies.  It is also safe to
# source from the real-data analysis: main() is called only when this file is
# executed with Rscript.

# -----------------------------------------------------------------------------
# General utilities
# -----------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

WORKFLOW_VERSION <- "2026-08-31-gaussian-pc5-v1"

RESULT_SCHEMA <- c(
  "study", "generator", "J", "n", "n_signal", "outcome", "relationship",
  "rho", "replicate", "method", "q", "fdp", "power", "n_selected"
)

q_key <- function(q) paste0("q", sprintf("%a", as.numeric(q)))

clip_probability <- function(x, eps = 1e-8) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

safe_r2 <- function(y, prediction) {
  ok <- is.finite(y) & is.finite(prediction)
  y <- as.numeric(y[ok])
  prediction <- as.numeric(prediction[ok])
  if (length(y) < 3L) return(NA_real_)
  denominator <- sum((y - mean(y))^2)
  if (!is.finite(denominator) || denominator <= .Machine$double.eps) {
    return(NA_real_)
  }
  1 - sum((y - prediction)^2) / denominator
}

binary_cross_entropy <- function(y, probability, eps = 1e-8) {
  probability <- clip_probability(probability, eps = eps)
  -mean(y * log(probability) + (1 - y) * log(1 - probability))
}

standardize_safe <- function(x, tol = 1e-12) {
  x <- as.numeric(x)
  sx <- stats::sd(x)
  if (!is.finite(sx) || sx <= tol) return(rep(0, length(x)))
  (x - mean(x)) / sx
}

as_integer_scalar <- function(x, name, minimum = 1L) {
  numeric_value <- suppressWarnings(as.numeric(x))
  if (length(numeric_value) != 1L || !is.finite(numeric_value) ||
      numeric_value != floor(numeric_value) || numeric_value > .Machine$integer.max ||
      numeric_value < minimum) {
    stop(name, " must be an integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(numeric_value)
}

as_flag <- function(x, name = "flag") {
  if (is.logical(x) && length(x) == 1L && !is.na(x)) return(x)
  value <- tolower(trimws(as.character(x)))
  if (value %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (value %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop(name, " must be true or false.", call. = FALSE)
}

split_csv <- function(x) {
  if (is.null(x) || length(x) == 0L || !nzchar(trimws(x))) {
    return(character(0))
  }
  trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE))
}

# A dependency-free, stable hash used only to derive deterministic RNG seeds
# and short checkpoint identifiers.  It is not intended for cryptography.
stable_hash <- function(...) {
  bytes <- utf8ToInt(paste(..., collapse = "|"))
  modulus <- 2147483629
  value <- 104729
  for (byte in bytes) value <- (value * 131 + byte) %% modulus
  as.numeric(value)
}

stable_seed <- function(base_seed = 1L, ...) {
  modulus <- 2147483629
  value <- (as.numeric(base_seed) + stable_hash(...)) %% modulus
  if (!is.finite(value) || value < 1) value <- 1
  as.integer(value)
}

current_script_path <- function() {
  definition_environment <- environment(current_script_path)
  if (exists(
    ".SIMULATION_SCRIPT_FILE",
    envir = definition_environment,
    inherits = FALSE
  )) {
    return(get(
      ".SIMULATION_SCRIPT_FILE",
      envir = definition_environment,
      inherits = FALSE
    ))
  }
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
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  }
  NA_character_
}

.SIMULATION_SCRIPT_FILE <- current_script_path()

find_project_root <- function(start = NULL) {
  if (is.null(start)) {
    script <- current_script_path()
    start <- if (is.na(script)) getwd() else dirname(script)
  }
  path <- normalizePath(start, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(path, "README.md")) || dir.exists(file.path(path, ".git"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) return(normalizePath(getwd(), mustWork = FALSE))
    path <- parent
  }
}

resolve_output_path <- function(path, project_root) {
  if (!nzchar(path)) stop("output path cannot be empty.", call. = FALSE)
  path <- path.expand(path)
  is_absolute <- grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
  normalizePath(
    if (is_absolute) path else file.path(project_root, path),
    mustWork = FALSE
  )
}

write_csv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(path), ".tmp_"), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(x, temporary, row.names = FALSE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) stop("Could not write ", path, call. = FALSE)
  invisible(path)
}

save_rds_atomic <- function(object, path, compress = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(path), ".tmp_"), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = compress)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) stop("Could not write ", path, call. = FALSE)
  invisible(path)
}

parallel_lapply_reproducible <- function(X, FUN, cores = 1L) {
  cores <- max(1L, as.integer(cores))
  if (.Platform$OS.type == "unix" && cores > 1L) {
    parallel::mclapply(
      X,
      FUN,
      mc.cores = min(cores, length(X)),
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
  } else {
    if (cores > 1L) {
      message("Forked parallelism is unavailable; using sequential execution.")
    }
    lapply(X, FUN)
  }
}

required_packages_for_methods <- function(methods) {
  packages <- character(0)
  if ("Delporte_MDS" %in% methods) packages <- c(packages, "xgboost")
  if ("Dai_MDS" %in% methods) packages <- c(packages, "glmnet")
  if (any(c("Knockoff", "KnockoffPlus") %in% methods)) {
    packages <- c(packages, "knockoff", "glmnet")
  }
  unique(packages)
}

workflow_code_hash <- function(path = .SIMULATION_SCRIPT_FILE) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("Cannot hash the simulation workflow because its script path is unknown.",
         call. = FALSE)
  }
  unname(tools::md5sum(path))
}

package_version_signature <- function(methods) {
  packages <- required_packages_for_methods(methods)
  versions <- vapply(
    packages,
    function(package) {
      if (!requireNamespace(package, quietly = TRUE)) return("missing")
      as.character(utils::packageVersion(package))
    },
    character(1L)
  )
  c(R = paste(R.version$major, R.version$minor, sep = "."), versions)
}

check_required_packages <- function(methods, stop_on_missing = TRUE) {
  packages <- required_packages_for_methods(methods)
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L && stop_on_missing) {
    stop(
      "Missing package(s) required by the selected methods: ",
      paste(missing, collapse = ", "),
      ". Install them before running; this script never installs packages.",
      call. = FALSE
    )
  }
  missing
}

has_method_error <- function(error) {
  length(error) == 1L && !is.na(error) && nzchar(trimws(as.character(error)))
}

summarize_split_errors <- function(diagnostics, method) {
  failed <- which(!is.na(diagnostics$error) & nzchar(trimws(diagnostics$error)))
  if (length(failed) == 0L) return(NA_character_)
  messages <- unique(diagnostics$error[failed])
  paste0(
    method, " failed on ", length(failed), " of ", nrow(diagnostics),
    " splits. First error: ", messages[[1L]]
  )
}

# -----------------------------------------------------------------------------
# Data generation
# -----------------------------------------------------------------------------

rdirichlet_base <- function(n, alpha) {
  n <- as.integer(n)
  alpha <- as.numeric(alpha)
  if (n < 1L || length(alpha) < 2L || any(!is.finite(alpha)) || any(alpha <= 0)) {
    stop("Invalid Dirichlet dimensions or concentration parameters.", call. = FALSE)
  }
  gamma_draws <- matrix(
    stats::rgamma(n * length(alpha), shape = rep(alpha, each = n), rate = 1),
    nrow = n,
    ncol = length(alpha)
  )
  gamma_draws <- pmax(gamma_draws, .Machine$double.xmin)
  gamma_draws / rowSums(gamma_draws)
}

generate_single_dirichlet <- function(n, J, alpha = 0.5, seed = 1L) {
  set.seed(seed)
  X <- rdirichlet_base(n, rep(alpha, J))
  colnames(X) <- paste0("X", seq_len(J))
  list(X = X, component = rep(1L, n), generator = "single_dirichlet")
}

generate_two_component_dirichlet <- function(
    n,
    J,
    alpha_high = 0.75,
    alpha_low = 0.25,
    mixture_probability = 0.5,
    seed = 1L) {

  if (J < 2L || mixture_probability <= 0 || mixture_probability >= 1) {
    stop("Invalid two-component Dirichlet configuration.", call. = FALSE)
  }
  set.seed(seed)
  split_point <- floor(J / 2L)
  alpha_1 <- c(rep(alpha_high, split_point), rep(alpha_low, J - split_point))
  alpha_2 <- c(rep(alpha_low, split_point), rep(alpha_high, J - split_point))
  component <- stats::rbinom(n, size = 1L, prob = mixture_probability) + 1L
  X <- matrix(NA_real_, nrow = n, ncol = J)
  first <- which(component == 1L)
  second <- which(component == 2L)
  if (length(first) > 0L) X[first, ] <- rdirichlet_base(length(first), alpha_1)
  if (length(second) > 0L) X[second, ] <- rdirichlet_base(length(second), alpha_2)
  colnames(X) <- paste0("X", seq_len(J))
  list(X = X, component = component, generator = "two_component_dirichlet")
}

generate_latent_block_gaussian <- function(
    n,
    J,
    rho,
    block_size = 10L,
    seed = 1L) {

  if (J %% block_size != 0L) stop("J must be divisible by block_size.", call. = FALSE)
  if (!is.finite(rho) || rho < 0 || rho >= 1) stop("rho must be in [0, 1).", call. = FALSE)
  set.seed(seed)
  Z <- matrix(NA_real_, nrow = n, ncol = J)
  n_blocks <- J %/% block_size
  for (block in seq_len(n_blocks)) {
    columns <- ((block - 1L) * block_size + 1L):(block * block_size)
    common_factor <- stats::rnorm(n)
    independent_noise <- matrix(stats::rnorm(n * block_size), nrow = n)
    Z[, columns] <- sqrt(rho) * common_factor + sqrt(1 - rho) * independent_noise
  }

  # Stable row-wise softmax.  The design parameter rho is the latent Gaussian
  # within-block correlation; closure changes correlations observed in X.
  row_maximum <- apply(Z, 1L, max)
  exponentiated <- exp(Z - row_maximum)
  X <- exponentiated / rowSums(exponentiated)
  colnames(X) <- paste0("X", seq_len(J))
  colnames(Z) <- paste0("Z", seq_len(J))
  list(
    X = X,
    Z_latent = Z,
    component = rep(1L, n),
    generator = "latent_block_gaussian_softmax",
    rho = rho,
    block_size = block_size
  )
}

pairwise_balance_features <- function(X, eps = 1e-12) {
  J <- ncol(X)
  if (J %% 2L != 0L) stop("J must be even for paired balances.", call. = FALSE)
  odd <- seq.int(1L, J, by = 2L)
  even <- odd + 1L
  balances <- (
    log(pmax(X[, odd, drop = FALSE], eps)) -
      log(pmax(X[, even, drop = FALSE], eps))
  ) / sqrt(2)
  colnames(balances) <- paste0("B", seq_len(ncol(balances)))
  balances
}

log_abundance <- function(X, eps = 1e-12) log(pmax(as.matrix(X), eps))

model_matrix_for_relationship <- function(X, relationship = c("linear", "nonlinear")) {
  relationship <- match.arg(relationship)
  if (relationship == "linear") log_abundance(X) else pairwise_balance_features(X)
}

perturbed_model_matrix_no_closure <- function(
    X,
    Z,
    j,
    replacement,
    relationship = c("linear", "nonlinear"),
    eps = 1e-12) {

  relationship <- match.arg(relationship)
  Z_new <- Z
  if (relationship == "linear") {
    Z_new[, j] <- log(pmax(replacement, eps))
    return(Z_new)
  }

  pair <- (j + 1L) %/% 2L
  if (j %% 2L == 1L) {
    partner <- j + 1L
    Z_new[, pair] <- (
      log(pmax(replacement, eps)) - log(pmax(X[, partner], eps))
    ) / sqrt(2)
  } else {
    partner <- j - 1L
    Z_new[, pair] <- (
      log(pmax(X[, partner], eps)) - log(pmax(replacement, eps))
    ) / sqrt(2)
  }
  Z_new
}

make_signal_pair_order <- function(J, seed = 1L) {
  if (J %% 2L != 0L) stop("J must be even.", call. = FALSE)
  set.seed(seed)
  sample.int(J %/% 2L, size = J %/% 2L, replace = FALSE)
}

signal_features_from_pairs <- function(signal_pairs) {
  sort(as.integer(c(rbind(2L * signal_pairs - 1L, 2L * signal_pairs))))
}

generate_outcome <- function(
    X,
    pair_order,
    n_signal,
    outcome = c("continuous", "binary"),
    relationship = c("linear", "nonlinear"),
    target_r2,
    seed = 1L) {

  outcome <- match.arg(outcome)
  relationship <- match.arg(relationship)
  if (n_signal < 2L || n_signal %% 2L != 0L) {
    stop("n_signal must be a positive even number.", call. = FALSE)
  }
  signal_pairs <- pair_order[seq_len(n_signal %/% 2L)]
  truth <- signal_features_from_pairs(signal_pairs)

  if (relationship == "linear") {
    beta <- numeric(ncol(X))
    for (pair in signal_pairs) {
      beta[2L * pair - 1L] <- 1
      beta[2L * pair] <- -1
    }
    eta <- drop(log_abundance(X) %*% beta)
  } else {
    balances <- pairwise_balance_features(X)
    terms <- balances[, signal_pairs, drop = FALSE]^2
    terms <- apply(terms, 2L, standardize_safe)
    if (is.null(dim(terms))) terms <- matrix(terms, ncol = 1L)
    eta <- rowSums(terms)
    beta <- rep(NA_real_, ncol(X))
  }

  eta <- eta - mean(eta)
  variance_eta <- stats::var(eta)
  if (!is.finite(variance_eta) || variance_eta <= 0 ||
      !is.finite(target_r2) || target_r2 <= 0 || target_r2 >= 1) {
    stop("Outcome signal variance and target_r2 must be positive and finite.", call. = FALSE)
  }
  noise_variance <- variance_eta * (1 - target_r2) / target_r2
  set.seed(seed)
  latent_y <- eta + stats::rnorm(nrow(X), sd = sqrt(noise_variance))

  if (outcome == "continuous") {
    y <- latent_y
    cutpoint <- NA_real_
  } else {
    cutpoint <- stats::median(latent_y)
    y <- as.integer(latent_y > cutpoint)
  }

  list(
    y = y,
    latent_y = latent_y,
    eta = eta,
    truth = truth,
    signal_pairs = signal_pairs,
    beta = beta,
    target_r2 = target_r2,
    oracle_latent_r2 = safe_r2(latent_y, eta),
    noise_variance = noise_variance,
    dichotomization_cutpoint = cutpoint,
    case_fraction = if (outcome == "binary") mean(y) else NA_real_
  )
}

# -----------------------------------------------------------------------------
# Gaussian raw-PC negative controls
# -----------------------------------------------------------------------------

make_balanced_fold_id <- function(n, nfolds, seed = 1L, strata = NULL) {
  set.seed(seed)
  fold_id <- integer(n)
  if (is.null(strata)) {
    fold_id[sample.int(n)] <- rep(seq_len(nfolds), length.out = n)
  } else {
    for (level in unique(strata)) {
      rows <- which(strata == level)
      fold_id[sample(rows, length(rows), replace = FALSE)] <-
        rep(seq_len(nfolds), length.out = length(rows))
    }
  }
  fold_id
}

generate_gaussian_pc_controls <- function(
    X,
    K = 5L,
    fold_id = NULL,
    seed = 1L,
    support = c("clip", "none"),
    lower = 0,
    upper = 1,
    verbose = FALSE) {

  support <- match.arg(support)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X)
  J <- ncol(X)
  K <- as_integer_scalar(K, "K")
  if (n < 6L || J < 2L || any(!is.finite(X))) {
    stop("X must be a finite numeric matrix with n >= 6 and J >= 2.", call. = FALSE)
  }
  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
    stop("lower and upper must be finite with lower < upper.", call. = FALSE)
  }
  if (is.null(fold_id)) fold_id <- make_balanced_fold_id(n, 2L, seed)
  fold_id <- as.integer(fold_id)
  if (length(fold_id) != n || !identical(sort(unique(fold_id)), 1:2)) {
    stop("fold_id must assign every row to fold 1 or 2.", call. = FALSE)
  }

  controls <- matrix(NA_real_, nrow = n, ncol = J)
  diagnostics <- vector("list", 2L * J)
  diagnostic_index <- 0L

  for (j in seq_len(J)) {
    remaining <- setdiff(seq_len(J), j)
    for (heldout_fold in 1:2) {
      train_rows <- which(fold_id != heldout_fold)
      test_rows <- which(fold_id == heldout_fold)
      raw_train <- X[train_rows, remaining, drop = FALSE]
      raw_test <- X[test_rows, remaining, drop = FALSE]
      K_effective <- min(K, nrow(raw_train) - 1L, ncol(raw_train))
      fallback <- FALSE

      fitted_values <- tryCatch({
        # Deliberately raw coordinates: no log, CLR, ILR, remainder
        # normalization, or compositional transformation is used here.
        pca <- stats::prcomp(
          raw_train,
          center = TRUE,
          scale. = FALSE,
          rank. = K_effective
        )
        train_scores <- pca$x[, seq_len(K_effective), drop = FALSE]
        test_scores <- stats::predict(pca, newdata = raw_test)[
          , seq_len(K_effective), drop = FALSE
        ]
        train_design <- cbind(`(Intercept)` = 1, train_scores)
        test_design <- cbind(`(Intercept)` = 1, test_scores)
        fit <- stats::lm.fit(train_design, X[train_rows, j])
        coefficients <- as.numeric(fit$coefficients)
        coefficients[!is.finite(coefficients)] <- 0
        prediction <- drop(test_design %*% coefficients)
        rank <- fit$rank %||% ncol(train_design)
        residual_sd <- sqrt(sum(fit$residuals^2) / max(length(train_rows) - rank, 1L))
        list(prediction = prediction, residual_sd = residual_sd)
      }, error = function(e) NULL)

      if (is.null(fitted_values) || any(!is.finite(fitted_values$prediction))) {
        fallback <- TRUE
        fitted_values <- list(
          prediction = rep(mean(X[train_rows, j]), length(test_rows)),
          residual_sd = stats::sd(X[train_rows, j])
        )
      }
      residual_sd <- fitted_values$residual_sd
      if (!is.finite(residual_sd) || residual_sd < 0) residual_sd <- 0

      set.seed(stable_seed(seed, "gaussian_pc", j, heldout_fold))
      raw_draw <- stats::rnorm(
        length(test_rows),
        mean = fitted_values$prediction,
        sd = residual_sd
      )
      below <- sum(raw_draw < lower)
      above <- sum(raw_draw > upper)

      # Clipping enforces the declared numerical abundance support only.  It
      # does not rescale any other feature and is therefore not closure.
      draw <- if (support == "clip") pmin(pmax(raw_draw, lower), upper) else raw_draw
      controls[test_rows, j] <- draw

      diagnostic_index <- diagnostic_index + 1L
      diagnostics[[diagnostic_index]] <- data.frame(
        feature = j,
        heldout_fold = heldout_fold,
        n_generated = length(test_rows),
        K_effective = K_effective,
        residual_sd = residual_sd,
        n_below_support = below,
        n_above_support = above,
        n_clipped = if (support == "clip") below + above else 0L,
        clip_rate = if (support == "clip") (below + above) / length(test_rows) else 0,
        n_at_or_below_log_floor = sum(draw <= 1e-12),
        used_fallback = fallback,
        stringsAsFactors = FALSE
      )
    }
    if (verbose && (j %% 25L == 0L || j == J)) {
      message("Generated Gaussian PC controls for ", j, " of ", J, " features.")
    }
  }

  if (any(!is.finite(controls))) stop("Non-finite negative controls were generated.", call. = FALSE)
  colnames(controls) <- colnames(X) %||% paste0("X", seq_len(J))
  diagnostics <- do.call(rbind, diagnostics)
  rownames(diagnostics) <- NULL
  list(
    # `controls` is the descriptive public name; `values` is retained as a
    # concise compatibility alias for callers that already consume it.
    controls = controls,
    values = controls,
    fold_id = fold_id,
    diagnostics = diagnostics,
    K = K,
    support = support,
    lower = lower,
    upper = upper,
    downstream_log_floor = 1e-12,
    transformation = "centered_raw_X_minus_j",
    predictive_model = "Gaussian linear model on first K PC scores",
    closure_adjustment = FALSE
  )
}

replace_feature_no_closure <- function(X, j, replacement) {
  X <- as.matrix(X)
  if (length(replacement) != nrow(X)) stop("replacement has the wrong length.", call. = FALSE)
  perturbed <- X
  perturbed[, j] <- replacement
  # No row normalization or adjustment to X_-j is performed.
  perturbed
}

# -----------------------------------------------------------------------------
# Mirror selection and MDS aggregation
# -----------------------------------------------------------------------------

mirror_statistic <- function(beta_1, beta_2) {
  statistic <- sign(beta_1 * beta_2) * (abs(beta_1) + abs(beta_2))
  statistic[!is.finite(statistic)] <- 0
  statistic
}

select_by_mirror <- function(statistic, q, offset = 1L) {
  offset <- as_integer_scalar(offset, "mirror offset", minimum = 0L)
  thresholds <- sort(unique(abs(statistic[is.finite(statistic) & statistic != 0])))
  if (length(thresholds) == 0L) return(integer(0))
  for (threshold in thresholds) {
    estimated_fdp <- (offset + sum(statistic <= -threshold, na.rm = TRUE)) /
      max(sum(statistic >= threshold, na.rm = TRUE), 1L)
    if (estimated_fdp <= q) return(which(statistic >= threshold))
  }
  integer(0)
}

aggregate_mds_selection <- function(selection_list, J, q) {
  inclusion_rate <- numeric(J)
  if (length(selection_list) == 0L) {
    return(list(selected = integer(0), inclusion_rate = inclusion_rate))
  }
  for (selected in selection_list) {
    selected <- unique(as.integer(selected))
    selected <- selected[selected >= 1L & selected <= J]
    if (length(selected) > 0L) inclusion_rate[selected] <-
      inclusion_rate[selected] + 1 / length(selected)
  }
  inclusion_rate <- inclusion_rate / length(selection_list)
  if (all(inclusion_rate <= 0)) {
    return(list(selected = integer(0), inclusion_rate = inclusion_rate))
  }
  ordering <- order(inclusion_rate, decreasing = FALSE)
  admissible <- which(cumsum(inclusion_rate[ordering]) <= q)
  threshold <- if (length(admissible) == 0L) 0 else
    inclusion_rate[ordering[max(admissible)]]
  list(selected = which(inclusion_rate > threshold), inclusion_rate = inclusion_rate)
}

selection_metrics <- function(selected, truth, J) {
  selected <- unique(as.integer(selected))
  selected <- selected[selected >= 1L & selected <= J]
  true_positive <- sum(selected %in% truth)
  false_positive <- length(selected) - true_positive
  list(
    fdp = false_positive / max(length(selected), 1L),
    power = true_positive / max(length(truth), 1L),
    n_selected = length(selected)
  )
}

empty_selections <- function(methods, q_values) {
  setNames(
    lapply(methods, function(method) {
      setNames(lapply(q_values, function(q) integer(0)), vapply(q_values, q_key, character(1)))
    }),
    methods
  )
}

# -----------------------------------------------------------------------------
# Proposed Delporte MDS (XGBoost)
# -----------------------------------------------------------------------------

default_xgboost_control <- function() {
  list(
    nrounds = 350L,
    max_depth = 2L,
    eta = 0.04,
    min_child_weight = 2,
    subsample = 0.85,
    colsample_bytree = 0.85,
    lambda = 2,
    alpha = 0,
    nthread = 1L
  )
}

fit_xgboost_outcome <- function(Z, y, outcome, control, seed) {
  if (outcome == "binary" && length(unique(y)) < 2L) {
    return(list(kind = "constant", value = mean(y), outcome = outcome))
  }
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(Z), label = y)
  parameters <- list(
    objective = if (outcome == "continuous") "reg:squarederror" else "binary:logistic",
    eval_metric = if (outcome == "continuous") "rmse" else "logloss",
    max_depth = control$max_depth,
    eta = control$eta,
    min_child_weight = control$min_child_weight,
    subsample = control$subsample,
    colsample_bytree = control$colsample_bytree,
    lambda = control$lambda,
    alpha = control$alpha,
    nthread = control$nthread,
    seed = as.integer(seed),
    verbosity = 0L
  )
  model <- xgboost::xgb.train(
    params = parameters,
    data = dtrain,
    nrounds = as.integer(control$nrounds),
    verbose = 0L
  )
  list(kind = "xgboost", model = model, outcome = outcome)
}

predict_xgboost_outcome <- function(fit, Z) {
  if (fit$kind == "constant") return(rep(fit$value, nrow(Z)))
  as.numeric(stats::predict(fit$model, newdata = as.matrix(Z)))
}

make_rotating_plan <- function(n, D = 50L, seed = 1L, strata = NULL) {
  D <- as_integer_scalar(D, "D")
  rounds <- ceiling(D / 5L)
  plan <- vector("list", rounds * 5L)
  index <- 0L
  for (round in seq_len(rounds)) {
    folds <- make_balanced_fold_id(n, 10L, stable_seed(seed, "round", round), strata)
    for (rotation in seq_len(5L)) {
      first_fold <- 2L * rotation - 1L
      second_fold <- 2L * rotation
      index <- index + 1L
      plan[[index]] <- list(
        train = which(!folds %in% c(first_fold, second_fold)),
        test1 = which(folds == first_fold),
        test2 = which(folds == second_fold),
        round = round,
        rotation = rotation
      )
    }
  }
  plan[seq_len(D)]
}

run_delporte_mds <- function(
    X,
    y,
    outcome = c("continuous", "binary"),
    relationship = c("linear", "nonlinear"),
    negative_controls = NULL,
    q_values = c(0.05, 0.10),
    D = 50L,
    xgb_control = default_xgboost_control(),
    K = 5L,
    support = "clip",
    seed = 1L,
    verbose = FALSE) {

  outcome <- match.arg(outcome)
  relationship <- match.arg(relationship)
  X <- as.matrix(X)
  y <- as.numeric(y)
  J <- ncol(X)
  if (is.null(negative_controls)) {
    negative_controls <- generate_gaussian_pc_controls(
      X, K = K, seed = stable_seed(seed, "negative_controls"), support = support
    )
  }
  control_object <- if (is.list(negative_controls)) negative_controls else
    list(values = as.matrix(negative_controls), diagnostics = NULL)
  T_tilde <- as.matrix(control_object$values)
  if (!identical(dim(T_tilde), dim(X))) {
    stop("negative_controls must have the same dimensions as X.", call. = FALSE)
  }

  split_selections <- setNames(
    lapply(q_values, function(q) vector("list", D)),
    vapply(q_values, q_key, character(1))
  )
  mirror_matrix <- matrix(0, nrow = D, ncol = J)
  beta_1_matrix <- matrix(0, nrow = D, ncol = J)
  beta_2_matrix <- matrix(0, nrow = D, ncol = J)
  diagnostics <- vector("list", D)
  strata <- if (outcome == "binary") y else NULL
  plan <- make_rotating_plan(nrow(X), D = D, seed = stable_seed(seed, "plan"), strata = strata)
  Z <- model_matrix_for_relationship(X, relationship)

  for (split in seq_len(D)) {
    split_result <- tryCatch({
      indices <- plan[[split]]
      fit <- fit_xgboost_outcome(
        Z[indices$train, , drop = FALSE], y[indices$train], outcome,
        xgb_control, stable_seed(seed, "xgboost", split)
      )
      prediction_1 <- predict_xgboost_outcome(fit, Z[indices$test1, , drop = FALSE])
      prediction_2 <- predict_xgboost_outcome(fit, Z[indices$test2, , drop = FALSE])
      baseline_1 <- if (outcome == "continuous") {
        safe_r2(y[indices$test1], prediction_1)
      } else {
        binary_cross_entropy(y[indices$test1], prediction_1)
      }
      baseline_2 <- if (outcome == "continuous") {
        safe_r2(y[indices$test2], prediction_2)
      } else {
        binary_cross_entropy(y[indices$test2], prediction_2)
      }

      beta_1 <- numeric(J)
      beta_2 <- numeric(J)
      X_1 <- X[indices$test1, , drop = FALSE]
      X_2 <- X[indices$test2, , drop = FALSE]
      Z_test1 <- Z[indices$test1, , drop = FALSE]
      Z_test2 <- Z[indices$test2, , drop = FALSE]
      for (j in seq_len(J)) {
        Z_1 <- perturbed_model_matrix_no_closure(
          X_1,
          Z_test1,
          j,
          T_tilde[indices$test1, j],
          relationship
        )
        Z_2 <- perturbed_model_matrix_no_closure(
          X_2,
          Z_test2,
          j,
          T_tilde[indices$test2, j],
          relationship
        )
        changed_prediction_1 <- predict_xgboost_outcome(fit, Z_1)
        changed_prediction_2 <- predict_xgboost_outcome(fit, Z_2)
        changed_1 <- if (outcome == "continuous") {
          safe_r2(y[indices$test1], changed_prediction_1)
        } else {
          binary_cross_entropy(y[indices$test1], changed_prediction_1)
        }
        changed_2 <- if (outcome == "continuous") {
          safe_r2(y[indices$test2], changed_prediction_2)
        } else {
          binary_cross_entropy(y[indices$test2], changed_prediction_2)
        }
        if (outcome == "continuous") {
          beta_1[j] <- baseline_1 - changed_1
          beta_2[j] <- baseline_2 - changed_2
        } else {
          beta_1[j] <- changed_1 - baseline_1
          beta_2[j] <- changed_2 - baseline_2
        }
      }
      statistic <- mirror_statistic(beta_1, beta_2)
      list(
        beta_1 = beta_1,
        beta_2 = beta_2,
        statistic = statistic,
        baseline_1 = baseline_1,
        baseline_2 = baseline_2,
        error = NA_character_
      )
    }, error = function(e) {
      list(
        beta_1 = numeric(J), beta_2 = numeric(J), statistic = numeric(J),
        baseline_1 = NA_real_, baseline_2 = NA_real_,
        error = conditionMessage(e)
      )
    })

    beta_1_matrix[split, ] <- split_result$beta_1
    beta_2_matrix[split, ] <- split_result$beta_2
    mirror_matrix[split, ] <- split_result$statistic
    for (q in q_values) {
      split_selections[[q_key(q)]][[split]] <-
        select_by_mirror(split_result$statistic, q, offset = 0L)
    }
    diagnostics[[split]] <- data.frame(
      split = split,
      baseline_test1 = split_result$baseline_1,
      baseline_test2 = split_result$baseline_2,
      n_nonzero_mirrors = sum(split_result$statistic != 0),
      error = split_result$error,
      stringsAsFactors = FALSE
    )
    if (verbose && (split %% 10L == 0L || split == D)) {
      message("Completed Delporte split ", split, " of ", D, ".")
    }
  }

  final <- setNames(vector("list", length(q_values)), vapply(q_values, q_key, character(1)))
  rates <- setNames(vector("list", length(q_values)), vapply(q_values, q_key, character(1)))
  for (q in q_values) {
    aggregated <- aggregate_mds_selection(split_selections[[q_key(q)]], J, q)
    final[[q_key(q)]] <- aggregated$selected
    rates[[q_key(q)]] <- aggregated$inclusion_rate
  }
  split_diagnostics <- do.call(rbind, diagnostics)
  list(
    selections = list(Delporte_MDS = final),
    inclusion_rates = list(Delporte_MDS = rates),
    split_selections = split_selections,
    split_diagnostics = split_diagnostics,
    mirror_diagnostics = list(
      beta_1 = beta_1_matrix,
      beta_2 = beta_2_matrix,
      statistic = mirror_matrix
    ),
    negative_control_diagnostics = control_object$diagnostics,
    learner = "XGBoost",
    perturbation = "GaussianPC5_NoClosure",
    mirror_offset = 0L,
    error = summarize_split_errors(split_diagnostics, "Delporte MDS")
  )
}

# -----------------------------------------------------------------------------
# Dai MDS comparator
# -----------------------------------------------------------------------------

fit_dai_first_half <- function(Z, y, outcome, nfolds, seed) {
  if (nrow(Z) < 3L) stop("Dai MDS requires at least three training rows.")
  nfolds <- max(3L, min(as.integer(nfolds), nrow(Z)))
  strata <- if (outcome == "binary") y else NULL
  fold_id <- make_balanced_fold_id(nrow(Z), nfolds, seed, strata)
  family <- if (outcome == "continuous") "gaussian" else "binomial"
  type_measure <- if (outcome == "continuous") "mse" else "deviance"
  cv_fit <- glmnet::cv.glmnet(
    x = Z,
    y = y,
    family = family,
    alpha = 1,
    nfolds = nfolds,
    foldid = fold_id,
    type.measure = type_measure,
    standardize = TRUE,
    intercept = TRUE,
    parallel = FALSE
  )
  coefficients <- as.numeric(stats::coef(cv_fit, s = "lambda.min"))[-1L]
  coefficients[!is.finite(coefficients)] <- 0
  coefficients
}

fit_dai_second_half <- function(Z, y, screen, J, outcome) {
  coefficients <- numeric(J)
  if (length(screen) == 0L) return(coefficients)
  if (outcome == "continuous") {
    fit <- stats::lm.fit(Z[, screen, drop = FALSE], y)
  } else {
    fit <- suppressWarnings(stats::glm.fit(
      x = Z[, screen, drop = FALSE],
      y = y,
      family = stats::binomial(),
      intercept = FALSE,
      singular.ok = TRUE
    ))
  }
  values <- as.numeric(fit$coefficients)
  values[!is.finite(values)] <- 0
  coefficients[screen] <- values
  coefficients
}

make_half_split <- function(n, seed, strata = NULL) {
  if (is.null(strata)) {
    set.seed(seed)
    first <- sample.int(n, size = floor(n / 2L), replace = FALSE)
  } else {
    first <- integer(0)
    for (level in unique(strata)) {
      rows <- which(strata == level)
      set.seed(stable_seed(seed, "stratum", level))
      first <- c(first, sample(rows, size = floor(length(rows) / 2L), replace = FALSE))
    }
  }
  list(first = sort(first), second = setdiff(seq_len(n), first))
}

run_dai_mds <- function(
    Z,
    y,
    outcome = c("continuous", "binary"),
    q_values = c(0.05, 0.10),
    D = 50L,
    nfolds = 10L,
    seed = 1L,
    verbose = FALSE) {

  outcome <- match.arg(outcome)
  Z <- as.matrix(Z)
  y <- as.numeric(y)
  J <- ncol(Z)
  split_selections <- setNames(
    lapply(q_values, function(q) vector("list", D)),
    vapply(q_values, q_key, character(1))
  )
  mirror_matrix <- matrix(0, nrow = D, ncol = J)
  diagnostics <- vector("list", D)

  for (split in seq_len(D)) {
    result <- tryCatch({
      split_rows <- make_half_split(
        nrow(Z), stable_seed(seed, "dai_split", split),
        if (outcome == "binary") y else NULL
      )
      beta_1 <- fit_dai_first_half(
        Z[split_rows$first, , drop = FALSE], y[split_rows$first], outcome,
        nfolds, stable_seed(seed, "dai_cv", split)
      )
      screen <- which(beta_1 != 0)
      beta_2 <- fit_dai_second_half(
        Z[split_rows$second, , drop = FALSE], y[split_rows$second],
        screen, J, outcome
      )
      list(statistic = mirror_statistic(beta_1, beta_2), n_screened = length(screen),
           error = NA_character_)
    }, error = function(e) {
      list(statistic = numeric(J), n_screened = 0L, error = conditionMessage(e))
    })
    mirror_matrix[split, ] <- result$statistic
    for (q in q_values) {
      split_selections[[q_key(q)]][[split]] <-
        select_by_mirror(result$statistic, q, offset = 0L)
    }
    diagnostics[[split]] <- data.frame(
      split = split,
      n_screened = result$n_screened,
      n_nonzero_mirrors = sum(result$statistic != 0),
      error = result$error,
      stringsAsFactors = FALSE
    )
    if (verbose && (split %% 10L == 0L || split == D)) {
      message("Completed Dai split ", split, " of ", D, ".")
    }
  }

  final <- setNames(vector("list", length(q_values)), vapply(q_values, q_key, character(1)))
  rates <- setNames(vector("list", length(q_values)), vapply(q_values, q_key, character(1)))
  for (q in q_values) {
    aggregated <- aggregate_mds_selection(split_selections[[q_key(q)]], J, q)
    final[[q_key(q)]] <- aggregated$selected
    rates[[q_key(q)]] <- aggregated$inclusion_rate
  }
  split_diagnostics <- do.call(rbind, diagnostics)
  list(
    selections = list(Dai_MDS = final),
    inclusion_rates = list(Dai_MDS = rates),
    split_selections = split_selections,
    split_diagnostics = split_diagnostics,
    mirror_diagnostics = list(statistic = mirror_matrix),
    mirror_offset = 0L,
    error = summarize_split_errors(split_diagnostics, "Dai MDS")
  )
}

# -----------------------------------------------------------------------------
# Univariate BY comparator
# -----------------------------------------------------------------------------

univariate_continuous_pvalues <- function(Z, y) {
  centered_Z <- sweep(Z, 2L, colMeans(Z), FUN = "-")
  centered_y <- y - mean(y)
  denominator <- sqrt(sum(centered_y^2) * colSums(centered_Z^2))
  correlation <- drop(crossprod(centered_y, centered_Z)) / denominator
  correlation[!is.finite(correlation)] <- 0
  correlation <- pmax(pmin(correlation, 1 - 1e-14), -1 + 1e-14)
  degrees_freedom <- nrow(Z) - 2L
  statistic <- correlation * sqrt(
    degrees_freedom / pmax(1 - correlation^2, .Machine$double.eps)
  )
  p_values <- 2 * stats::pt(-abs(statistic), df = degrees_freedom)
  p_values[!is.finite(p_values)] <- 1
  p_values
}

univariate_binary_pvalues <- function(Z, y) {
  group_0 <- which(y == 0)
  group_1 <- which(y == 1)
  if (length(group_0) < 2L || length(group_1) < 2L) return(rep(1, ncol(Z)))
  Z_0 <- Z[group_0, , drop = FALSE]
  Z_1 <- Z[group_1, , drop = FALSE]
  mean_0 <- colMeans(Z_0)
  mean_1 <- colMeans(Z_1)
  variance_0 <- colSums(sweep(Z_0, 2L, mean_0, FUN = "-")^2) / (length(group_0) - 1L)
  variance_1 <- colSums(sweep(Z_1, 2L, mean_1, FUN = "-")^2) / (length(group_1) - 1L)
  component_0 <- variance_0 / length(group_0)
  component_1 <- variance_1 / length(group_1)
  standard_error_squared <- component_0 + component_1
  statistic <- (mean_1 - mean_0) / sqrt(standard_error_squared)
  degrees_freedom <- standard_error_squared^2 / (
    component_0^2 / (length(group_0) - 1L) +
      component_1^2 / (length(group_1) - 1L)
  )
  p_values <- 2 * stats::pt(-abs(statistic), df = degrees_freedom)
  p_values[!is.finite(p_values)] <- 1
  pmin(pmax(p_values, 0), 1)
}

run_univariate_by <- function(
    Z,
    y,
    outcome = c("continuous", "binary"),
    q_values = c(0.05, 0.10)) {

  outcome <- match.arg(outcome)
  Z <- as.matrix(Z)
  y <- as.numeric(y)
  failure <- NA_character_
  p_values <- tryCatch(
    if (outcome == "continuous") univariate_continuous_pvalues(Z, y) else
      univariate_binary_pvalues(Z, y),
    error = function(e) {
      failure <<- conditionMessage(e)
      rep(NA_real_, ncol(Z))
    }
  )
  adjusted <- if (has_method_error(failure)) {
    rep(NA_real_, length(p_values))
  } else {
    stats::p.adjust(p_values, method = "BY")
  }
  selected <- setNames(
    lapply(q_values, function(q) {
      if (has_method_error(failure)) integer(0) else which(adjusted <= q)
    }),
    vapply(q_values, q_key, character(1))
  )
  list(
    selections = list(Univariate_BY = selected),
    p_values = p_values,
    adjusted_p_values = adjusted,
    test = if (outcome == "continuous") "univariate linear regression" else
      "Welch two-sample t-test",
    error = failure
  )
}

# -----------------------------------------------------------------------------
# Knockoff and Knockoff+ comparators
# -----------------------------------------------------------------------------

deterministic_glmnet_coefdiff <- function(
    X,
    Xk,
    y,
    outcome,
    nfolds,
    seed,
    nlambda = 500L) {

  X <- as.matrix(X)
  Xk <- as.matrix(Xk)
  if (!identical(dim(X), dim(Xk))) {
    stop("X and Xk must have identical dimensions.", call. = FALSE)
  }
  p <- ncol(X)
  set.seed(stable_seed(seed, "swap"))
  swap <- stats::rbinom(p, size = 1L, prob = 0.5)
  swap_matrix <- matrix(swap, nrow = nrow(X), ncol = p, byrow = TRUE)
  X_swapped <- X * (1 - swap_matrix) + Xk * swap_matrix
  Xk_swapped <- X * swap_matrix + Xk * (1 - swap_matrix)

  design <- cbind(X_swapped, Xk_swapped)
  design <- sweep(design, 2L, colMeans(design), FUN = "-")
  scales <- apply(design, 2L, stats::sd)
  scales[!is.finite(scales) | scales <= .Machine$double.eps] <- 1
  design <- sweep(design, 2L, scales, FUN = "/")

  family <- if (outcome == "continuous") "gaussian" else "binomial"
  strata <- if (outcome == "binary") y else NULL
  if (nrow(design) < 3L) stop("Knockoff requires at least three rows.")
  nfolds <- max(3L, min(as.integer(nfolds), nrow(design)))
  fold_id <- make_balanced_fold_id(
    nrow(design),
    nfolds,
    seed = stable_seed(seed, "folds"),
    strata = strata
  )
  fit_arguments <- list(
    x = design,
    y = y,
    family = family,
    intercept = TRUE,
    standardize = FALSE,
    standardize.response = FALSE,
    parallel = FALSE,
    foldid = fold_id,
    nfolds = nfolds
  )
  if (outcome == "continuous") {
    lambda_max <- max(abs(drop(crossprod(design, y)))) / nrow(design)
    if (!is.finite(lambda_max) || lambda_max <= .Machine$double.eps) {
      return(numeric(p))
    }
    nlambda <- as_integer_scalar(nlambda, "nlambda")
    exponent <- (seq_len(nlambda) - 1L) / nlambda
    fit_arguments$lambda <- lambda_max * (1 / 2000)^exponent
  }

  fit <- do.call(glmnet::cv.glmnet, fit_arguments)
  coefficients <- as.numeric(stats::coef(fit, s = "lambda.min"))[-1L]
  if (length(coefficients) != 2L * p) {
    stop("glmnet returned an unexpected coefficient vector.", call. = FALSE)
  }
  original <- seq_len(p)
  statistic <- abs(coefficients[original]) - abs(coefficients[original + p])
  statistic * (1 - 2 * swap)
}

run_knockoffs <- function(
    Z,
    y,
    outcome = c("continuous", "binary"),
    q_values = c(0.05, 0.10),
    seed = 1L,
    Xk = NULL,
    nfolds = 10L) {

  outcome <- match.arg(outcome)
  Z <- as.matrix(Z)
  y <- as.numeric(y)
  methods <- c("Knockoff", "KnockoffPlus")
  failure <- NA_character_
  threshold_errors <- character(0)
  statistic <- tryCatch({
    set.seed(stable_seed(seed, "knockoff_copy"))
    if (is.null(Xk)) {
      Xk <- knockoff::create.second_order(Z, method = "equi", shrink = TRUE)
    }
    deterministic_glmnet_coefdiff(
      X = Z,
      Xk = Xk,
      y = y,
      outcome = outcome,
      nfolds = nfolds,
      seed = stable_seed(seed, "knockoff_statistic")
    )
  }, error = function(e) {
    failure <<- conditionMessage(e)
    numeric(ncol(Z))
  })

  selections <- empty_selections(methods, q_values)
  if (any(statistic != 0)) {
    for (q in q_values) {
      threshold_0 <- tryCatch(
        knockoff::knockoff.threshold(statistic, fdr = q, offset = 0),
        error = function(e) {
          threshold_errors <<- c(
            threshold_errors,
            paste0("Knockoff threshold at q=", q, ": ", conditionMessage(e))
          )
          Inf
        }
      )
      threshold_1 <- tryCatch(
        knockoff::knockoff.threshold(statistic, fdr = q, offset = 1),
        error = function(e) {
          threshold_errors <<- c(
            threshold_errors,
            paste0("Knockoff+ threshold at q=", q, ": ", conditionMessage(e))
          )
          Inf
        }
      )
      selections$Knockoff[[q_key(q)]] <-
        if (is.finite(threshold_0)) which(statistic >= threshold_0) else integer(0)
      selections$KnockoffPlus[[q_key(q)]] <-
        if (is.finite(threshold_1)) which(statistic >= threshold_1) else integer(0)
    }
  }
  if (!has_method_error(failure) && length(threshold_errors) > 0L) {
    failure <- paste(unique(threshold_errors), collapse = " | ")
  }
  list(
    selections = selections,
    W = statistic,
    Xk = Xk,
    error = failure,
    offsets = c(Knockoff = 0L, KnockoffPlus = 1L)
  )
}

# -----------------------------------------------------------------------------
# Scenario definitions
# -----------------------------------------------------------------------------

main_target_r2 <- function(J, n_signal, outcome) {
  if (J == 100L) return(if (outcome == "continuous") 0.20 else 0.30)
  if (J == 1000L && outcome == "continuous") {
    return(if (n_signal == 10L) 0.30 else 0.60)
  }
  if (J == 1000L && outcome == "binary") {
    return(if (n_signal == 10L) 0.45 else 0.90)
  }
  stop("No target R2 is defined for this scenario.", call. = FALSE)
}

make_main_scenarios <- function(
    J = c(100L, 1000L),
    outcome = c("continuous", "binary"),
    generator = "single_dirichlet") {

  rows <- list()
  index <- 0L
  for (dimension in as.integer(J)) {
    signal_counts <- if (dimension == 100L) c(6L, 10L) else c(10L, 20L)
    for (n_signal in signal_counts) {
      for (outcome_type in outcome) {
        for (relationship in c("linear", "nonlinear")) {
          index <- index + 1L
          rows[[index]] <- data.frame(
            scenario_id = paste(
              "main", generator, paste0("J", dimension), paste0("s", n_signal),
              outcome_type, relationship, sep = "_"
            ),
            scenario_number = index,
            study = "main",
            generator = generator,
            J = dimension,
            n = 1000L,
            n_signal = n_signal,
            outcome = outcome_type,
            relationship = relationship,
            rho = NA_real_,
            target_r2 = main_target_r2(dimension, n_signal, outcome_type),
            dirichlet_alpha = if (generator == "single_dirichlet") 0.5 else NA_real_,
            mixture_alpha_high = if (generator == "two_component_dirichlet") 0.75 else NA_real_,
            mixture_alpha_low = if (generator == "two_component_dirichlet") 0.25 else NA_real_,
            mixture_probability = if (generator == "two_component_dirichlet") 0.5 else NA_real_,
            block_size = NA_integer_,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  scenarios <- do.call(rbind, rows)
  scenarios$scenario_number <- seq_len(nrow(scenarios))
  scenarios
}

make_correlation_scenarios <- function() {
  grid <- expand.grid(
    n_signal = c(6L, 10L),
    relationship = c("linear", "nonlinear"),
    rho = c(0, 0.3, 0.6),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$scenario_number <- seq_len(nrow(grid))
  grid$scenario_id <- with(
    grid,
    paste0(
      "correlation_latent_block_gaussian_softmax_J100_s", n_signal,
      "_continuous_", relationship, "_rho", gsub("\\.", "", format(rho))
    )
  )
  grid$study <- "correlation"
  grid$generator <- "latent_block_gaussian_softmax"
  grid$J <- 100L
  grid$n <- 1000L
  grid$outcome <- "continuous"
  grid$target_r2 <- 0.20
  grid$dirichlet_alpha <- NA_real_
  grid$mixture_alpha_high <- NA_real_
  grid$mixture_alpha_low <- NA_real_
  grid$mixture_probability <- NA_real_
  grid$block_size <- 10L
  grid[, c(
    "scenario_id", "scenario_number", "study", "generator", "J", "n",
    "n_signal", "outcome", "relationship", "rho", "target_r2",
    "dirichlet_alpha", "mixture_alpha_high", "mixture_alpha_low",
    "mixture_probability", "block_size"
  )]
}

validate_scenarios <- function(scenarios) {
  required <- c(
    "scenario_id", "scenario_number", "study", "generator", "J", "n",
    "n_signal", "outcome", "relationship", "rho", "target_r2"
  )
  missing <- setdiff(required, names(scenarios))
  if (length(missing) > 0L) {
    stop("Scenario table is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(scenarios$scenario_id)) stop("scenario_id values must be unique.", call. = FALSE)
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# One replicate and result normalization
# -----------------------------------------------------------------------------

canonical_method_names <- function(x) {
  if (length(x) == 0L || any(tolower(x) == "all")) {
    return(c("Delporte_MDS", "Dai_MDS", "Univariate_BY", "Knockoff", "KnockoffPlus"))
  }
  aliases <- c(
    delporte = "Delporte_MDS",
    delporte_mds = "Delporte_MDS",
    dai = "Dai_MDS",
    dai_mds = "Dai_MDS",
    by = "Univariate_BY",
    univariate_by = "Univariate_BY",
    knockoff = "Knockoff",
    knockoff_plus = "KnockoffPlus",
    knockoffplus = "KnockoffPlus"
  )
  keys <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  unknown <- keys[!keys %in% names(aliases)]
  if (length(unknown) > 0L) {
    stop("Unknown method(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  unique(unname(aliases[keys]))
}

composition_group_key <- function(scenario) {
  paste(
    scenario$study,
    scenario$generator,
    scenario$J,
    ifelse(is.na(scenario$rho), "NA", scenario$rho),
    sep = "|"
  )
}

generate_composition_for_scenario <- function(scenario, replicate, seed_base) {
  seed <- stable_seed(
    seed_base, "composition", scenario$study, scenario$generator,
    scenario$J, scenario$rho, replicate
  )
  switch(
    scenario$generator,
    single_dirichlet = generate_single_dirichlet(
      scenario$n, scenario$J, alpha = scenario$dirichlet_alpha, seed = seed
    ),
    two_component_dirichlet = generate_two_component_dirichlet(
      scenario$n,
      scenario$J,
      alpha_high = scenario$mixture_alpha_high,
      alpha_low = scenario$mixture_alpha_low,
      mixture_probability = scenario$mixture_probability,
      seed = seed
    ),
    latent_block_gaussian_softmax = generate_latent_block_gaussian(
      scenario$n, scenario$J, scenario$rho,
      block_size = scenario$block_size, seed = seed
    ),
    stop("Unknown generator: ", scenario$generator, call. = FALSE)
  )
}

normalize_selection_rows <- function(
    selections,
    requested_methods,
    q_values,
    truth,
    scenario,
    replicate,
    failed_methods = character(0)) {

  rows <- list()
  index <- 0L
  for (method in requested_methods) {
    method_selections <- selections[[method]] %||%
      setNames(lapply(q_values, function(q) integer(0)), vapply(q_values, q_key, character(1)))
    for (q in q_values) {
      metrics <- if (method %in% failed_methods) {
        list(fdp = NA_real_, power = NA_real_, n_selected = NA_integer_)
      } else {
        selection_metrics(
          method_selections[[q_key(q)]] %||% integer(0),
          truth,
          scenario$J
        )
      }
      index <- index + 1L
      rows[[index]] <- data.frame(
        study = as.character(scenario$study),
        generator = as.character(scenario$generator),
        J = as.integer(scenario$J),
        n = as.integer(scenario$n),
        n_signal = as.integer(scenario$n_signal),
        outcome = as.character(scenario$outcome),
        relationship = as.character(scenario$relationship),
        rho = as.numeric(scenario$rho),
        replicate = as.integer(replicate),
        method = method,
        q = as.numeric(q),
        fdp = metrics$fdp,
        power = metrics$power,
        n_selected = metrics$n_selected,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

run_one_replicate <- function(replicate, scenarios, config) {
  groups <- split(seq_len(nrow(scenarios)), vapply(
    seq_len(nrow(scenarios)),
    function(i) composition_group_key(scenarios[i, , drop = FALSE]),
    character(1)
  ))
  metric_rows <- list()
  control_diagnostics <- list()
  method_diagnostics <- list()
  row_index <- 0L
  diagnostic_index <- 0L
  method_diagnostic_index <- 0L

  for (group_rows in groups) {
    representative <- scenarios[group_rows[1L], , drop = FALSE]
    composition <- generate_composition_for_scenario(representative, replicate, config$seed_base)
    X <- composition$X
    group_seed <- stable_seed(
      config$seed_base, "group", composition_group_key(representative), replicate
    )
    pair_order <- make_signal_pair_order(ncol(X), stable_seed(group_seed, "pairs"))

    negative_controls <- NULL
    if ("Delporte_MDS" %in% config$methods) {
      negative_controls <- generate_gaussian_pc_controls(
        X,
        K = config$K,
        seed = stable_seed(group_seed, "controls"),
        support = config$support,
        verbose = FALSE
      )
      diagnostic <- negative_controls$diagnostics
      diagnostic$study <- representative$study
      diagnostic$generator <- representative$generator
      diagnostic$J <- representative$J
      diagnostic$rho <- representative$rho
      diagnostic$replicate <- replicate
      diagnostic_index <- diagnostic_index + 1L
      control_diagnostics[[diagnostic_index]] <- diagnostic
    }

    Z_comparator <- log_abundance(X)
    Xk <- NULL
    if (any(c("Knockoff", "KnockoffPlus") %in% config$methods)) {
      Xk <- tryCatch({
        set.seed(stable_seed(group_seed, "shared_knockoff_copy"))
        knockoff::create.second_order(Z_comparator, method = "equi", shrink = TRUE)
      }, error = function(e) NULL)
    }

    for (scenario_row in group_rows) {
      scenario <- scenarios[scenario_row, , drop = FALSE]
      scenario_seed <- stable_seed(
        config$seed_base, "scenario", scenario$scenario_id, replicate
      )
      response <- generate_outcome(
        X = X,
        pair_order = pair_order,
        n_signal = scenario$n_signal,
        outcome = scenario$outcome,
        relationship = scenario$relationship,
        target_r2 = scenario$target_r2,
        seed = stable_seed(scenario_seed, "outcome")
      )
      selections <- list()
      method_errors <- setNames(rep(NA_character_, length(config$methods)), config$methods)

      if ("Delporte_MDS" %in% config$methods) {
        result <- run_delporte_mds(
          X = X,
          y = response$y,
          outcome = scenario$outcome,
          relationship = scenario$relationship,
          negative_controls = negative_controls,
          q_values = config$q_values,
          D = config$n_splits,
          xgb_control = config$xgb_control,
          K = config$K,
          support = config$support,
          seed = stable_seed(scenario_seed, "delporte")
        )
        selections$Delporte_MDS <- result$selections$Delporte_MDS
        if (has_method_error(result$error)) {
          method_errors[["Delporte_MDS"]] <- result$error
        }
      }

      if ("Dai_MDS" %in% config$methods) {
        result <- run_dai_mds(
          Z = Z_comparator,
          y = response$y,
          outcome = scenario$outcome,
          q_values = config$q_values,
          D = config$n_splits,
          nfolds = config$dai_nfolds,
          seed = stable_seed(scenario_seed, "dai")
        )
        selections$Dai_MDS <- result$selections$Dai_MDS
        if (has_method_error(result$error)) {
          method_errors[["Dai_MDS"]] <- result$error
        }
      }

      if ("Univariate_BY" %in% config$methods) {
        result <- run_univariate_by(
          Z = Z_comparator,
          y = response$y,
          outcome = scenario$outcome,
          q_values = config$q_values
        )
        selections$Univariate_BY <- result$selections$Univariate_BY
        if (has_method_error(result$error)) {
          method_errors[["Univariate_BY"]] <- result$error
        }
      }

      if (any(c("Knockoff", "KnockoffPlus") %in% config$methods)) {
        result <- run_knockoffs(
          Z = Z_comparator,
          y = response$y,
          outcome = scenario$outcome,
          q_values = config$q_values,
          seed = stable_seed(scenario_seed, "knockoff"),
          Xk = Xk,
          nfolds = config$knockoff_nfolds
        )
        if ("Knockoff" %in% config$methods) {
          selections$Knockoff <- result$selections$Knockoff
        }
        if ("KnockoffPlus" %in% config$methods) {
          selections$KnockoffPlus <- result$selections$KnockoffPlus
        }
        if (has_method_error(result$error)) {
          affected <- intersect(c("Knockoff", "KnockoffPlus"), config$methods)
          method_errors[affected] <- result$error
        }
      }

      missing_results <- setdiff(config$methods, names(selections))
      if (length(missing_results) > 0L) {
        method_errors[missing_results] <- "Method returned no selection object."
      }
      failed_methods <- names(method_errors)[
        vapply(method_errors, has_method_error, logical(1L))
      ]
      method_diagnostic_index <- method_diagnostic_index + 1L
      method_diagnostics[[method_diagnostic_index]] <- data.frame(
        study = as.character(scenario$study),
        scenario_id = as.character(scenario$scenario_id),
        generator = as.character(scenario$generator),
        J = as.integer(scenario$J),
        n_signal = as.integer(scenario$n_signal),
        outcome = as.character(scenario$outcome),
        relationship = as.character(scenario$relationship),
        rho = as.numeric(scenario$rho),
        replicate = as.integer(replicate),
        method = config$methods,
        status = ifelse(config$methods %in% failed_methods, "failed", "success"),
        error = unname(method_errors[config$methods]),
        stringsAsFactors = FALSE
      )

      row_index <- row_index + 1L
      metric_rows[[row_index]] <- normalize_selection_rows(
        selections = selections,
        requested_methods = config$methods,
        q_values = config$q_values,
        truth = response$truth,
        scenario = scenario,
        replicate = replicate,
        failed_methods = failed_methods
      )
    }
  }

  list(
    replicate = replicate,
    metrics = do.call(rbind, metric_rows),
    method_diagnostics = do.call(rbind, method_diagnostics),
    negative_control_diagnostics = if (length(control_diagnostics) > 0L) {
      do.call(rbind, control_diagnostics)
    } else {
      data.frame()
    }
  )
}

checkpoint_path <- function(config, replicate) {
  file.path(
    config$output_dir,
    "checkpoints",
    config$run_id,
    sprintf("replicate_%03d.rds", replicate)
  )
}

run_one_replicate_checkpointed <- function(replicate, scenarios, config) {
  path <- checkpoint_path(config, replicate)
  if (file.exists(path)) {
    existing <- tryCatch(readRDS(path), error = function(e) NULL)
    valid <- !is.null(existing) &&
      identical(existing$run_signature, config$run_signature) &&
      identical(existing$replicate, as.integer(replicate)) &&
      is.data.frame(existing$metrics) &&
      identical(names(existing$metrics), RESULT_SCHEMA) &&
      is.data.frame(existing$method_diagnostics) &&
      nrow(existing$method_diagnostics) > 0L &&
      all(existing$method_diagnostics$status == "success")
    if (isTRUE(valid)) {
      return(existing)
    }
  }
  result <- run_one_replicate(replicate, scenarios, config)
  result$run_signature <- config$run_signature
  result$completed_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  save_rds_atomic(result, path, compress = TRUE)
  result
}

summarize_replicate_metrics <- function(metrics) {
  group_columns <- c(
    "study", "generator", "J", "n", "n_signal", "outcome",
    "relationship", "rho", "method", "q"
  )
  key_parts <- lapply(metrics[group_columns], function(x) {
    ifelse(is.na(x), "NA", as.character(x))
  })
  group_key <- do.call(paste, c(key_parts, sep = "|"))
  rows <- lapply(split(metrics, group_key), function(d) {
    fdp <- d$fdp[is.finite(d$fdp)]
    power <- d$power[is.finite(d$power)]
    selected <- d$n_selected[is.finite(d$n_selected)]
    safe_mean <- function(x) if (length(x) > 0L) mean(x) else NA_real_
    safe_median <- function(x) if (length(x) > 0L) stats::median(x) else NA_real_
    safe_quantile <- function(x, probability) {
      if (length(x) > 0L) {
        as.numeric(stats::quantile(x, probability, names = FALSE))
      } else {
        NA_real_
      }
    }
    safe_mcse <- function(x) {
      if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
    }
    data.frame(
      study = d$study[1L],
      generator = d$generator[1L],
      J = d$J[1L],
      n = d$n[1L],
      n_signal = d$n_signal[1L],
      outcome = d$outcome[1L],
      relationship = d$relationship[1L],
      rho = d$rho[1L],
      method = d$method[1L],
      q = d$q[1L],
      n_replicates = nrow(d),
      n_success = length(fdp),
      n_failed = nrow(d) - length(fdp),
      mean_fdr = safe_mean(fdp),
      fdr_mcse = safe_mcse(fdp),
      median_fdp = safe_median(fdp),
      fdp_q025 = safe_quantile(fdp, 0.025),
      fdp_q975 = safe_quantile(fdp, 0.975),
      mean_power = safe_mean(power),
      power_mcse = safe_mcse(power),
      median_power = safe_median(power),
      power_q025 = safe_quantile(power, 0.025),
      power_q975 = safe_quantile(power, 0.975),
      mean_selected = safe_mean(selected),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  summary
}

# -----------------------------------------------------------------------------
# CLI and study runner
# -----------------------------------------------------------------------------

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  parsed <- list()
  index <- 1L
  while (index <= length(args)) {
    token <- args[[index]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token, call. = FALSE)
    token <- substring(token, 3L)
    if (grepl("=", token, fixed = TRUE)) {
      pieces <- strsplit(token, "=", fixed = TRUE)[[1L]]
      key <- pieces[[1L]]
      value <- paste(pieces[-1L], collapse = "=")
    } else {
      key <- token
      if (key %in% c("dry-run", "help")) {
        value <- "true"
      } else {
        if (index == length(args) || startsWith(args[[index + 1L]], "--")) {
          stop("Missing value for --", key, call. = FALSE)
        }
        index <- index + 1L
        value <- args[[index]]
      }
    }
    key <- gsub("-", "_", tolower(key), fixed = TRUE)
    parsed[[key]] <- value
    index <- index + 1L
  }
  parsed
}

environment_default <- function(primary, secondary = NULL, fallback = "") {
  value <- Sys.getenv(primary, unset = "")
  if (!nzchar(value) && !is.null(secondary)) value <- Sys.getenv(secondary, unset = "")
  if (!nzchar(value)) fallback else value
}

default_cores <- function() {
  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (length(detected) != 1L || !is.finite(detected)) return(1L)
  max(1L, min(10L, as.integer(detected) - 1L))
}

build_run_config <- function(cli = list(), project_root = find_project_root()) {
  supported_options <- c(
    "study", "j", "outcome", "replicates", "splits", "cores", "methods",
    "dry_run", "generator", "support", "output", "q", "k", "dai_nfolds",
    "knockoff_nfolds", "seed", "xgb_nrounds", "help"
  )
  unknown_options <- setdiff(names(cli), supported_options)
  if (length(unknown_options) > 0L) {
    stop(
      "Unknown option(s): --",
      paste(gsub("_", "-", unknown_options, fixed = TRUE), collapse = ", --"),
      ". Use --help to list supported options.",
      call. = FALSE
    )
  }
  study <- cli$study %||% environment_default("SIM_STUDY", fallback = "main")
  study <- match.arg(tolower(study), c("main", "correlation"))
  J_text <- cli$j %||% environment_default("SIM_J", fallback = "")
  outcome_text <- cli$outcome %||% environment_default("SIM_OUTCOME", fallback = "")
  replicates <- cli$replicates %||%
    environment_default("SIM_REPLICATES", "N_REPLICATES", "200")
  splits <- cli$splits %||% environment_default("SIM_SPLITS", "N_SPLITS", "50")
  cores <- cli$cores %||% environment_default("SIM_CORES", "CORES", as.character(default_cores()))
  methods_text <- cli$methods %||% environment_default("SIM_METHODS", fallback = "all")
  dry_run_text <- cli$dry_run %||% environment_default("SIM_DRY_RUN", fallback = "false")
  generator_default <- if (study == "main") {
    "single_dirichlet"
  } else {
    "latent_block_gaussian_softmax"
  }
  generator <- cli$generator %||%
    environment_default("SIM_GENERATOR", fallback = generator_default)
  support <- cli$support %||% environment_default("SIM_NC_SUPPORT", fallback = "clip")
  output_text <- cli$output %||% environment_default("SIM_OUTPUT", "OUTPUT_DIR", "")

  all_dimensions <- !nzchar(J_text) || identical(tolower(trimws(J_text)), "all")
  J <- if (all_dimensions) {
    if (study == "main") c(100L, 1000L) else 100L
  } else {
    J_numeric <- suppressWarnings(as.numeric(split_csv(J_text)))
    if (any(!is.finite(J_numeric)) || any(J_numeric != floor(J_numeric))) {
      stop("--J must contain integer dimensions 100 and/or 1000.", call. = FALSE)
    }
    as.integer(J_numeric)
  }
  if (anyNA(J) || any(!J %in% c(100L, 1000L))) {
    stop("--J must contain 100 and/or 1000.", call. = FALSE)
  }
  all_outcomes <- !nzchar(outcome_text) ||
    identical(tolower(trimws(outcome_text)), "all")
  outcomes <- if (all_outcomes) {
    if (study == "main") c("continuous", "binary") else "continuous"
  } else {
    tolower(split_csv(outcome_text))
  }
  if (any(!outcomes %in% c("continuous", "binary"))) {
    stop("--outcome must be continuous and/or binary.", call. = FALSE)
  }
  if (study == "correlation" && (!identical(J, 100L) || !identical(outcomes, "continuous"))) {
    stop("The correlation study is defined only for J=100 continuous outcomes.", call. = FALSE)
  }
  generator <- if (study == "main") {
    match.arg(generator, c("single_dirichlet", "two_component_dirichlet"))
  } else {
    match.arg(generator, "latent_block_gaussian_softmax")
  }
  support <- match.arg(support, c("clip", "none"))
  methods <- canonical_method_names(split_csv(methods_text))
  q_values <- if (!is.null(cli$q)) as.numeric(split_csv(cli$q)) else c(0.05, 0.10)
  if (length(q_values) < 1L || any(!is.finite(q_values)) ||
      any(q_values <= 0) || any(q_values >= 1)) {
    stop("All q values must be in (0, 1).", call. = FALSE)
  }

  config <- list(
    project_root = project_root,
    study = study,
    J = unique(J),
    outcomes = unique(outcomes),
    generator = generator,
    n_replicates = as_integer_scalar(replicates, "replicates"),
    n_splits = as_integer_scalar(splits, "splits"),
    cores = as_integer_scalar(cores, "cores"),
    output_dir = NA_character_,
    methods = methods,
    dry_run = as_flag(dry_run_text, "dry-run"),
    q_values = sort(unique(q_values)),
    K = as_integer_scalar(cli$k %||% 5L, "K"),
    support = support,
    dai_nfolds = as_integer_scalar(cli$dai_nfolds %||% 10L, "dai-nfolds", 3L),
    knockoff_nfolds = as_integer_scalar(
      cli$knockoff_nfolds %||% 10L,
      "knockoff-nfolds",
      3L
    ),
    seed_base = as_integer_scalar(cli$seed %||% 20260831L, "seed"),
    xgb_control = default_xgboost_control(),
    workflow_version = WORKFLOW_VERSION,
    workflow_code_hash = workflow_code_hash(),
    package_versions = package_version_signature(methods)
  )
  if (!is.null(cli$xgb_nrounds)) {
    config$xgb_control$nrounds <- as_integer_scalar(cli$xgb_nrounds, "xgb-nrounds")
  }
  signature_fields <- c(
    config$study, config$J, config$outcomes, config$generator,
    config$n_replicates, config$n_splits, config$methods,
    config$q_values, config$K, config$support, config$seed_base,
    unlist(config$xgb_control, use.names = TRUE),
    config$dai_nfolds, config$knockoff_nfolds,
    config$workflow_version, config$workflow_code_hash,
    paste(names(config$package_versions), config$package_versions, sep = "=")
  )
  config$run_signature <- paste(signature_fields, collapse = "|")
  config$run_id <- sprintf("run_%010.0f", stable_hash(config$run_signature))
  config$output_is_explicit <- nzchar(output_text)
  config$output_dir <- if (config$output_is_explicit) {
    resolve_output_path(output_text, project_root)
  } else {
    file.path(
      project_root,
      "reproducibility",
      "results",
      config$study,
      config$generator,
      config$run_id
    )
  }
  config
}

make_configured_scenarios <- function(config) {
  scenarios <- if (config$study == "main") {
    make_main_scenarios(
      J = config$J,
      outcome = config$outcomes,
      generator = config$generator
    )
  } else {
    make_correlation_scenarios()
  }
  validate_scenarios(scenarios)
  scenarios
}

print_dry_run <- function(config, scenarios) {
  missing <- check_required_packages(config$methods, stop_on_missing = FALSE)
  message("Dry run: no files or simulations will be created.")
  message("Study: ", config$study)
  message("Generator: ", config$generator)
  message("Scenarios: ", nrow(scenarios))
  message("Replicates: ", config$n_replicates, "; MDS splits: ", config$n_splits)
  message("Methods: ", paste(config$methods, collapse = ", "))
  message("Negative controls: raw PC", config$K, " Gaussian, support=", config$support,
          ", no closure")
  message("Output: ", config$output_dir)
  if (length(missing) > 0L) {
    message("Missing packages for an actual run: ", paste(missing, collapse = ", "))
  }
  print(scenarios, row.names = FALSE)
  invisible(list(config = config, scenarios = scenarios, missing_packages = missing))
}

initialize_output_directory <- function(config) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  configuration_path <- file.path(config$output_dir, "run_configuration.rds")
  if (file.exists(configuration_path)) {
    existing <- tryCatch(readRDS(configuration_path), error = function(e) NULL)
    if (is.null(existing) ||
        !identical(existing$run_signature, config$run_signature)) {
      stop(
        "The output directory belongs to a different or unreadable run: ",
        config$output_dir,
        ". Choose a new --output path.",
        call. = FALSE
      )
    }
    return(invisible(configuration_path))
  }

  known_outputs <- c(
    "scenario_definitions.csv", "replicate_metrics.csv", "summary.csv",
    "method_diagnostics.csv", "negative_control_diagnostics.csv",
    "run_manifest.rds", "checkpoints"
  )
  conflicts <- known_outputs[
    file.exists(file.path(config$output_dir, known_outputs)) |
      dir.exists(file.path(config$output_dir, known_outputs))
  ]
  if (length(conflicts) > 0L) {
    stop(
      "The output directory contains results without a compatible run lock: ",
      paste(conflicts, collapse = ", "),
      ". Choose a new --output path.",
      call. = FALSE
    )
  }
  save_rds_atomic(
    list(
      run_signature = config$run_signature,
      run_id = config$run_id,
      workflow_version = config$workflow_version,
      workflow_code_hash = config$workflow_code_hash,
      configuration = config,
      created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ),
    configuration_path
  )
  invisible(configuration_path)
}

run_simulation_study <- function(config) {
  scenarios <- make_configured_scenarios(config)
  if (isTRUE(config$dry_run)) return(print_dry_run(config, scenarios))
  check_required_packages(config$methods)
  initialize_output_directory(config)
  write_csv_atomic(scenarios, file.path(config$output_dir, "scenario_definitions.csv"))

  message(
    "Running ", config$study, " study: ", nrow(scenarios), " scenarios x ",
    config$n_replicates, " replicates; methods = ", paste(config$methods, collapse = ", ")
  )
  results <- parallel_lapply_reproducible(
    seq_len(config$n_replicates),
    function(replicate) {
      run_one_replicate_checkpointed(replicate, scenarios, config)
    },
    cores = config$cores
  )
  bad <- vapply(results, function(x) !is.list(x) || is.null(x$metrics), logical(1))
  if (any(bad)) stop("At least one replicate failed before producing a checkpoint.", call. = FALSE)

  metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
  rownames(metrics) <- NULL
  metrics <- metrics[, RESULT_SCHEMA, drop = FALSE]
  metrics <- metrics[order(
    metrics$study, metrics$J, metrics$n_signal, metrics$outcome,
    metrics$relationship, metrics$rho, metrics$replicate, metrics$method, metrics$q,
    na.last = TRUE
  ), , drop = FALSE]
  summary <- summarize_replicate_metrics(metrics)
  diagnostics <- do.call(rbind, lapply(results, function(x) x$negative_control_diagnostics))
  method_diagnostics <- do.call(rbind, lapply(results, `[[`, "method_diagnostics"))
  rownames(method_diagnostics) <- NULL

  write_csv_atomic(metrics, file.path(config$output_dir, "replicate_metrics.csv"))
  write_csv_atomic(summary, file.path(config$output_dir, "summary.csv"))
  write_csv_atomic(
    method_diagnostics,
    file.path(config$output_dir, "method_diagnostics.csv")
  )
  if (!is.null(diagnostics) && nrow(diagnostics) > 0L) {
    write_csv_atomic(
      diagnostics,
      file.path(config$output_dir, "negative_control_diagnostics.csv")
    )
  }
  manifest <- list(
    scenarios = scenarios,
    configuration = config,
    unified_result_schema = RESULT_SCHEMA,
    checkpoint_directory = file.path("checkpoints", config$run_id),
    output_files = c(
      scenario_definitions = "scenario_definitions.csv",
      run_configuration = "run_configuration.rds",
      replicate_metrics = "replicate_metrics.csv",
      summary = "summary.csv",
      method_diagnostics = "method_diagnostics.csv",
      negative_control_diagnostics = if (nrow(diagnostics) > 0L) {
        "negative_control_diagnostics.csv"
      } else {
        NA_character_
      }
    ),
    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    session_info = utils::sessionInfo()
  )
  save_rds_atomic(manifest, file.path(config$output_dir, "run_manifest.rds"))
  n_failures <- sum(method_diagnostics$status == "failed", na.rm = TRUE)
  if (n_failures > 0L) {
    warning(
      n_failures,
      " scenario-replicate method run(s) failed. Their metrics are NA; see ",
      file.path(config$output_dir, "method_diagnostics.csv"),
      call. = FALSE
    )
  }
  message("Completed. Results are in: ", config$output_dir)
  invisible(list(metrics = metrics, summary = summary, manifest = manifest))
}

print_usage <- function() {
  cat(paste(
    "Usage: Rscript reproducibility/simulation/run_simulations.R [options]",
    "",
    "Options (both --key=value and --key value are accepted):",
    "  --study main|correlation",
    "  --J all|100|1000|100,1000",
    "  --outcome all|continuous|binary|continuous,binary",
    "  --replicates INTEGER",
    "  --splits INTEGER",
    "  --cores INTEGER",
    "  --output PATH",
    "  --methods all|delporte,dai,by,knockoff,knockoff_plus",
    "  --generator single_dirichlet|two_component_dirichlet (main study)",
    "              latent_block_gaussian_softmax (correlation study; fixed)",
    "  --q VALUES (comma-separated; default 0.05,0.10)",
    "  --k INTEGER (raw PCs for Gaussian controls; default 5)",
    "  --seed INTEGER",
    "  --dai-nfolds INTEGER",
    "  --knockoff-nfolds INTEGER",
    "  --xgb-nrounds INTEGER",
    "  --support clip|none",
    "  --dry-run",
    "  --help",
    "",
    "Environment equivalents use the SIM_ prefix, for example SIM_STUDY,",
    "SIM_J, SIM_OUTCOME, SIM_REPLICATES, SIM_SPLITS, SIM_CORES, SIM_OUTPUT,",
    "SIM_METHODS, and SIM_DRY_RUN.",
    sep = "\n"
  ), "\n")
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cli <- parse_cli_args(args)
  if (!is.null(cli$help) && as_flag(cli$help, "help")) {
    print_usage()
    return(invisible(NULL))
  }
  config <- build_run_config(cli, project_root = find_project_root())
  run_simulation_study(config)
}

if (sys.nframe() == 0L) {
  main()
}
