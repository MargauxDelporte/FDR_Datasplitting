#!/usr/bin/env Rscript

# Focused verification for the current reproducibility workflow.
# This deliberately uses tiny dimensions and does not run a manuscript study.

smoke_script_path <- function() {
  file_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  if (length(file_argument) != 1L) {
    stop("Run this smoke test with Rscript.", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
}

smoke_file <- smoke_script_path()
repository_root <- normalizePath(
  file.path(dirname(smoke_file), "..", ".."),
  mustWork = TRUE
)
simulation_file <- file.path(
  repository_root,
  "reproducibility",
  "simulation",
  "run_simulations.R"
)
plot_file <- file.path(
  repository_root,
  "reproducibility",
  "simulation",
  "plot_simulation_results.R"
)
qin_file <- file.path(
  repository_root,
  "reproducibility",
  "real_data",
  "qin_2014_analysis.R"
)

source(simulation_file, local = FALSE)
invisible(check_required_packages(
  c("Delporte_MDS", "Dai_MDS", "Knockoff", "KnockoffPlus")
))

message("1/5 Checking data generators ...")
dirichlet <- generate_single_dirichlet(n = 80L, J = 12L, seed = 101L)
mixture <- generate_two_component_dirichlet(n = 80L, J = 12L, seed = 102L)
correlated <- generate_latent_block_gaussian(
  n = 80L,
  J = 12L,
  rho = 0.3,
  block_size = 3L,
  seed = 103L
)
stopifnot(
  identical(dim(dirichlet$X), c(80L, 12L)),
  identical(dim(mixture$X), c(80L, 12L)),
  identical(dim(correlated$X), c(80L, 12L)),
  max(abs(rowSums(dirichlet$X) - 1)) < 1e-10,
  max(abs(rowSums(mixture$X) - 1)) < 1e-10,
  max(abs(rowSums(correlated$X) - 1)) < 1e-10
)
mixture_scenarios <- make_main_scenarios(
  J = 100L,
  outcome = "continuous",
  generator = "two_component_dirichlet"
)
stopifnot(
  all(is.na(mixture_scenarios$dirichlet_alpha)),
  all(mixture_scenarios$mixture_alpha_high == 0.75),
  all(mixture_scenarios$mixture_alpha_low == 0.25),
  all(mixture_scenarios$mixture_probability == 0.5),
  length(unique(q_key(c(0.0501, 0.0502)))) == 2L
)

message("2/5 Checking Gaussian raw-PC controls and no-closure replacement ...")
controls <- generate_gaussian_pc_controls(
  dirichlet$X,
  K = 5L,
  seed = 201L,
  support = "clip"
)
stopifnot(
  identical(dim(controls$values), dim(dirichlet$X)),
  all(is.finite(controls$values)),
  all(controls$values >= 0 & controls$values <= 1),
  identical(controls$transformation, "centered_raw_X_minus_j"),
  identical(controls$closure_adjustment, FALSE),
  all(controls$diagnostics$K_effective == 5L),
  "n_at_or_below_log_floor" %in% names(controls$diagnostics)
)
perturbed <- replace_feature_no_closure(
  dirichlet$X,
  j = 1L,
  replacement = controls$values[, 1L]
)
stopifnot(
  identical(perturbed[, -1L, drop = FALSE], dirichlet$X[, -1L, drop = FALSE]),
  identical(perturbed[, 1L], controls$values[, 1L]),
  any(abs(rowSums(perturbed) - 1) > 1e-10)
)
Z_linear <- model_matrix_for_relationship(dirichlet$X, "linear")
Z_nonlinear <- model_matrix_for_relationship(dirichlet$X, "nonlinear")
stopifnot(
  isTRUE(all.equal(
    perturbed_model_matrix_no_closure(
      dirichlet$X, Z_linear, 1L, controls$values[, 1L], "linear"
    ),
    model_matrix_for_relationship(perturbed, "linear")
  )),
  isTRUE(all.equal(
    perturbed_model_matrix_no_closure(
      dirichlet$X, Z_nonlinear, 1L, controls$values[, 1L], "nonlinear"
    ),
    model_matrix_for_relationship(perturbed, "nonlinear")
  ))
)

message("3/5 Checking continuous and binary proposed-method branches ...")
pair_order <- make_signal_pair_order(12L, seed = 301L)
continuous <- generate_outcome(
  dirichlet$X,
  pair_order = pair_order,
  n_signal = 4L,
  outcome = "continuous",
  relationship = "linear",
  target_r2 = 0.30,
  seed = 302L
)
binary <- generate_outcome(
  dirichlet$X,
  pair_order = pair_order,
  n_signal = 4L,
  outcome = "binary",
  relationship = "nonlinear",
  target_r2 = 0.30,
  seed = 303L
)
xgb_control <- default_xgboost_control()
xgb_control$nrounds <- 2L
xgb_control$max_depth <- 1L

delporte_continuous <- run_delporte_mds(
  X = dirichlet$X,
  y = continuous$y,
  outcome = "continuous",
  relationship = "linear",
  negative_controls = controls,
  q_values = c(0.05, 0.10),
  D = 2L,
  xgb_control = xgb_control,
  seed = 304L
)
delporte_binary <- run_delporte_mds(
  X = dirichlet$X,
  y = binary$y,
  outcome = "binary",
  relationship = "nonlinear",
  negative_controls = controls,
  q_values = 0.10,
  D = 1L,
  xgb_control = xgb_control,
  seed = 305L
)
invalid_xgb_control <- xgb_control
invalid_xgb_control$max_depth <- -1L
delporte_failure <- run_delporte_mds(
  X = dirichlet$X,
  y = continuous$y,
  outcome = "continuous",
  relationship = "linear",
  negative_controls = controls,
  q_values = 0.05,
  D = 1L,
  xgb_control = invalid_xgb_control,
  seed = 306L
)
stopifnot(
  identical(names(delporte_continuous$selections), "Delporte_MDS"),
  nrow(delporte_continuous$split_diagnostics) == 2L,
  nrow(delporte_binary$split_diagnostics) == 1L,
  is.na(delporte_continuous$error),
  is.na(delporte_binary$error),
  has_method_error(delporte_failure$error),
  !is.na(delporte_failure$split_diagnostics$error[[1L]]),
  all(is.finite(binary$y)),
  identical(sort(unique(binary$y)), c(0L, 1L))
)

message("4/5 Checking co-located competing methods ...")
Z <- log_abundance(dirichlet$X)
dai <- run_dai_mds(
  Z = Z,
  y = continuous$y,
  outcome = "continuous",
  q_values = c(0.05, 0.10),
  D = 2L,
  nfolds = 3L,
  seed = 401L
)
by <- run_univariate_by(
  Z = Z,
  y = continuous$y,
  outcome = "continuous",
  q_values = c(0.05, 0.10)
)
knockoffs <- run_knockoffs(
  Z = Z,
  y = continuous$y,
  outcome = "continuous",
  q_values = c(0.05, 0.10),
  nfolds = 3L,
  seed = 402L
)
dai_binary <- run_dai_mds(
  Z = Z,
  y = binary$y,
  outcome = "binary",
  q_values = 0.10,
  D = 1L,
  nfolds = 3L,
  seed = 403L
)
by_binary <- run_univariate_by(
  Z = Z,
  y = binary$y,
  outcome = "binary",
  q_values = 0.10
)
knockoffs_binary <- run_knockoffs(
  Z = Z,
  y = binary$y,
  outcome = "binary",
  q_values = 0.10,
  nfolds = 3L,
  seed = 404L
)
stopifnot(
  identical(names(dai$selections), "Dai_MDS"),
  identical(names(by$selections), "Univariate_BY"),
  identical(names(knockoffs$selections), c("Knockoff", "KnockoffPlus")),
  length(by$p_values) == 12L,
  length(knockoffs$W) == 12L,
  is.na(dai$error),
  is.na(by$error),
  is.na(knockoffs$error),
  is.na(dai_binary$error),
  is.na(by_binary$error),
  is.na(knockoffs_binary$error),
  length(select_by_mirror(c(5, 0, 0), q = 0.10, offset = 1L)) == 0L,
  identical(select_by_mirror(c(5, 0, 0), q = 0.10, offset = 0L), 1L)
)

integrated_scenario <- data.frame(
  scenario_id = "smoke_single_dirichlet_J12_s4_continuous_linear",
  scenario_number = 1L,
  study = "main",
  generator = "single_dirichlet",
  J = 12L,
  n = 80L,
  n_signal = 4L,
  outcome = "continuous",
  relationship = "linear",
  rho = NA_real_,
  target_r2 = 0.30,
  dirichlet_alpha = 0.5,
  mixture_alpha_high = NA_real_,
  mixture_alpha_low = NA_real_,
  mixture_probability = NA_real_,
  block_size = NA_integer_,
  stringsAsFactors = FALSE
)
integrated_result <- run_one_replicate(
  replicate = 1L,
  scenarios = integrated_scenario,
  config = list(
    seed_base = 405L,
    methods = c(
      "Delporte_MDS", "Dai_MDS", "Univariate_BY", "Knockoff", "KnockoffPlus"
    ),
    K = 5L,
    support = "clip",
    q_values = 0.10,
    n_splits = 1L,
    xgb_control = xgb_control,
    dai_nfolds = 3L,
    knockoff_nfolds = 3L
  )
)
stopifnot(
  nrow(integrated_result$metrics) == 5L,
  nrow(integrated_result$method_diagnostics) == 5L,
  all(integrated_result$method_diagnostics$status == "success"),
  all(is.finite(integrated_result$metrics$fdp)),
  all(is.finite(integrated_result$metrics$power))
)

message("5/5 Checking the unified result schema and plotting ...")
scenario_main <- data.frame(
  study = "main",
  generator = "single_dirichlet",
  J = 100L,
  n = 1000L,
  n_signal = 6L,
  outcome = "continuous",
  relationship = "linear",
  rho = NA_real_,
  stringsAsFactors = FALSE
)
scenario_correlation <- scenario_main
scenario_correlation$study <- "correlation"
scenario_correlation$generator <- "latent_block_gaussian_softmax"
scenario_correlation$rho <- 0.3

main_rows <- normalize_selection_rows(
  selections = c(
    delporte_continuous$selections,
    dai$selections,
    by$selections,
    knockoffs$selections
  ),
  requested_methods = c(
    "Delporte_MDS", "Dai_MDS", "Univariate_BY", "Knockoff", "KnockoffPlus"
  ),
  q_values = c(0.05, 0.10),
  truth = continuous$truth,
  scenario = scenario_main,
  replicate = 1L
)
correlation_rows <- main_rows
correlation_rows$study <- scenario_correlation$study
correlation_rows$generator <- scenario_correlation$generator
correlation_rows$rho <- scenario_correlation$rho
metrics <- rbind(main_rows, correlation_rows)
expected_schema <- c(
  "study", "generator", "J", "n", "n_signal", "outcome", "relationship",
  "rho", "replicate", "method", "q", "fdp", "power", "n_selected"
)
stopifnot(identical(names(metrics), expected_schema))
failed_rows <- normalize_selection_rows(
  selections = by$selections,
  requested_methods = "Univariate_BY",
  q_values = 0.05,
  truth = continuous$truth,
  scenario = scenario_main,
  replicate = 2L,
  failed_methods = "Univariate_BY"
)
failed_summary <- summarize_replicate_metrics(failed_rows)
stopifnot(
  all(is.na(failed_rows[, c("fdp", "power", "n_selected")])),
  failed_summary$n_success == 0L,
  failed_summary$n_failed == 1L,
  is.na(failed_summary$mean_fdr)
)

correlation_config <- build_run_config(
  list(study = "correlation", methods = "by", dry_run = "true"),
  project_root = repository_root
)
stopifnot(
  identical(correlation_config$generator, "latent_block_gaussian_softmax"),
  identical(correlation_config$workflow_version, WORKFLOW_VERSION),
  is.character(correlation_config$workflow_code_hash),
  nchar(correlation_config$workflow_code_hash) == 32L,
  inherits(
    try(build_run_config(list(replicate = "1"), repository_root), silent = TRUE),
    "try-error"
  ),
  inherits(
    try(build_run_config(list(q = ""), repository_root), silent = TRUE),
    "try-error"
  ),
  inherits(try(as_integer_scalar("1.5", "test"), silent = TRUE), "try-error")
)

plot_environment <- new.env(parent = globalenv())
sys.source(plot_file, envir = plot_environment)
qin_environment <- new.env(parent = globalenv())
sys.source(qin_file, envir = qin_environment)
stopifnot(
  identical(
    normalizePath(plot_environment$.plot_script_file, mustWork = TRUE),
    normalizePath(plot_file, mustWork = TRUE)
  ),
  identical(
    normalizePath(qin_environment$QIN_SCRIPT_FILE, mustWork = TRUE),
    normalizePath(qin_file, mustWork = TRUE)
  )
)
temporary_root <- tempfile("fdr_smoke_")
dir.create(temporary_root, recursive = TRUE)
metrics_file <- file.path(temporary_root, "replicate_metrics.csv")
figure_directory <- file.path(temporary_root, "figures")
utils::write.csv(metrics, metrics_file, row.names = FALSE)
plot_result <- plot_environment$plot_simulation_results(
  input = metrics_file,
  output = figure_directory,
  repository_root = repository_root
)
failed_metrics_file <- file.path(temporary_root, "failed_metrics.csv")
utils::write.csv(rbind(metrics, failed_rows), failed_metrics_file, row.names = FALSE)
failed_plot <- try(
  plot_environment$plot_simulation_results(
    input = failed_metrics_file,
    output = file.path(temporary_root, "failed_figures"),
    repository_root = repository_root
  ),
  silent = TRUE
)
stopifnot(
  file.exists(file.path(figure_directory, "plot_summary.csv")),
  file.exists(file.path(figure_directory, "main_fdp.png")),
  file.exists(file.path(figure_directory, "main_power.pdf")),
  file.exists(file.path(figure_directory, "correlation_fdp.png")),
  file.exists(file.path(figure_directory, "correlation_power.pdf")),
  length(plot_result$files) == 9L,
  inherits(failed_plot, "try-error")
)

message("All focused smoke checks passed. No full simulation was run.")
