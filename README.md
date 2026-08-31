# FDR control by permutation and data splitting

This repository contains the code accompanying:

> Delporte, M., & Shi, Y. (2026). *False Discovery Rate Control via
> Permutation and Data Splitting in Complex Datasets: Application to Liver
> Function Microbiome Data.* Manuscript submitted for publication.

The current, portable implementation is under [`reproducibility/`](reproducibility/).
It contains the simulation studies used in the revision, plotting code, and the
Qin et al. real-data analysis. The older `Scenario/`, `Functions/`,
`Functions Dai/`, and `Case study/` directories are retained unchanged for
historical provenance; they are not sourced by the current workflow.

## Analyses in the current workflow

All simulations use `n = 1000`, 200 Monte Carlo replicates by default, 50 data
splits, and target FDR levels 0.05 and 0.10.

| Study | J | Signals | Outcome | Signal strength | Data generator |
|---|---:|---:|---|---|---|
| Main | 100 | 6, 10 | Continuous | oracle R-squared 0.20 | Dirichlet |
| Main | 100 | 6, 10 | Binary | latent R-squared 0.30 | Dirichlet |
| Main | 1000 | 10, 20 | Continuous | oracle R-squared 0.30, 0.60 | Dirichlet |
| Main | 1000 | 10, 20 | Binary | latent R-squared 0.45, 0.90 | Dirichlet |
| Correlation | 100 | 6, 10 | Continuous | oracle R-squared 0.20 | latent block Gaussian |

Both linear and nonlinear response mechanisms are included. The correlation
study uses exchangeable 10-feature latent Gaussian blocks with within-block
correlation 0, 0.3, or 0.6, followed by a row-wise softmax. Functions for both
single and two-component Dirichlet generation are included in the simulation
runner.

The unified runner places the proposed procedure and all competing procedures
in the same file. It evaluates:

- Delporte MDS with XGBoost, rotating 80/10/10 splits, and the revision study's
  `c = 0` mirror-FDP estimate;
- Dai MDS with Lasso screening, a no-intercept OLS or logistic refit, and
  `c = 0`;
- univariate testing with Benjamini-Yekutieli adjustment;
- second-order Knockoff and Knockoff+.

## Negative-control model

The current implementation uses the requested Gaussian predictive model. It is
response-free and two-fold cross-fitted. For every feature `j`:

1. PCA is fit to the centered, raw abundance columns `X[, -j]` in the training
   fold. No log-ratio, CLR, or ILR transformation is used.
2. `X[, j]` is regressed on the first five PC scores with a Gaussian linear
   model.
3. A held-out negative-control value is drawn from a normal distribution with
   the predicted conditional mean and fitted residual standard deviation.
4. Only column `j` is replaced. The remaining columns are unchanged and rows
   are not re-normalized, so there is no closure adjustment.

By default, Gaussian draws are clipped to the abundance support `[0, 1]` before
they enter models that use log abundance. This is a component-wise support
constraint only; it does not impose closure. The runner records the clipping
diagnostics and exposes the support policy as an argument. With
`--support none`, the draw itself is left unbounded, but a downstream
log-abundance model still floors values at `1e-12`; diagnostics count values
affected by that floor.

## Repository map

```text
reproducibility/
  simulation/
    run_simulations.R            simulation, proposed method, and comparators
    plot_simulation_results.R    figures and plot-level summaries
  real_data/
    qin_2014_analysis.R          Qin et al. albumin and cirrhosis analyses
  tests/
    smoke_test.R                 focused checks; no full simulation rerun
  results/                       generated locally and ignored by Git
```

Every entry point derives the repository root from its own file location. It
can therefore be called from any working directory and contains no user-specific
absolute paths.

## Software

The simulation and comparison methods require R and these CRAN packages:

```r
install.packages(c("glmnet", "xgboost", "knockoff", "ggplot2"))
```

The Qin analysis additionally requires Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("curatedMetagenomicData", "SummarizedExperiment"))
```

Scripts check required packages and stop with an informative message. They do
not install or update packages automatically.

## Focused verification

The smoke test uses tiny dimensions and only a few splits. It checks both data
generators, the raw-PC Gaussian controls, direct no-closure replacement, the
result schema, transparent method-failure handling, CLI validation, and
plotting without rerunning the full Monte Carlo study:

```sh
Rscript reproducibility/tests/smoke_test.R
```

To inspect the complete scenario plan and resolved configuration without
fitting any models:

```sh
Rscript reproducibility/simulation/run_simulations.R --study main --dry-run
Rscript reproducibility/simulation/run_simulations.R --study correlation --dry-run
```

## Run the simulations

Run all main J=100 and J=1000 continuous and binary Dirichlet scenarios:

```sh
Rscript reproducibility/simulation/run_simulations.R \
  --study main --J all --outcome all \
  --replicates 200 --splits 50 --cores 10 \
  --output reproducibility/results/main/single_dirichlet
```

Run a smaller subset, for example the J=100 binary scenarios:

```sh
Rscript reproducibility/simulation/run_simulations.R \
  --study main --J 100 --outcome binary \
  --replicates 200 --splits 50 --cores 10 \
  --output reproducibility/results/main/J100_binary
```

The default main generator is a single Dirichlet distribution with
`alpha = 0.5`. To run the included symmetric two-component mixture instead:

```sh
Rscript reproducibility/simulation/run_simulations.R \
  --study main --generator two_component_dirichlet \
  --J all --outcome all --replicates 200 --splits 50 --cores 10 \
  --output reproducibility/results/main/two_component_dirichlet
```

Run the J=100 latent-block Gaussian correlation study:

```sh
Rscript reproducibility/simulation/run_simulations.R \
  --study correlation --replicates 200 --splits 50 --cores 10 \
  --output reproducibility/results/correlation
```

Without `--output`, the runner creates a configuration-specific directory under
`reproducibility/results/<study>/<generator>/<run_id>/`. Completed replicate
checkpoints are reused, so an interrupted command can be run again safely.
Random-number seeds are determined by the study, scenario, and replicate
identifiers; changing the number of workers does not change the generated
datasets. Checkpoints are keyed to the full configuration, workflow version,
exact runner contents, R version, and relevant package versions. An explicit
output directory is locked to one configuration, preventing a later command
from silently overwriting a different study. Checkpoints containing a failed
method are retained for diagnosis but are not treated as reusable successes.

## Plot simulation results

```sh
Rscript reproducibility/simulation/plot_simulation_results.R \
  --input reproducibility/results/main/single_dirichlet/replicate_metrics.csv \
  --output reproducibility/results/main/single_dirichlet/figures

Rscript reproducibility/simulation/plot_simulation_results.R \
  --input reproducibility/results/correlation/replicate_metrics.csv \
  --output reproducibility/results/correlation/figures
```

The plotting script reads the unified replicate-level schema, writes
`plot_summary.csv`, and produces PDF and PNG figures for FDP and power. It can
also plot a user-supplied compatible `replicate_metrics.csv`. By default it
stops if method failures are present, avoiding selective-failure bias;
`--allow-failures` is an explicit opt-in and records the retained counts.

## Qin et al. real-data analysis

The analysis downloads `QinN_2014.relative_abundance` through
`curatedMetagenomicData` with `counts = FALSE`. It uses taxa only, applies a 5%
prevalence filter, converts percentages to proportions when necessary, and does
not re-close the matrix after filtering. Outcomes are serum albumin
(continuous) and cirrhosis versus healthy status (binary).

```sh
Rscript reproducibility/real_data/qin_2014_analysis.R \
  --splits 50 --output reproducibility/results/qin_2014
```

The output includes selected taxa, selection counts, MDS inclusion rates,
univariate p-values, split diagnostics, `method_errors.csv`,
`method_status.csv`, a complete RDS result object, and a run manifest with
session information. The manifest also records the exact dated
`curatedMetagenomicData` resource selected at runtime, an input fingerprint,
and both analysis-workflow code identifiers.

## Output schema

Simulation replicate results use one long table:

```text
study, generator, J, n, n_signal, outcome, relationship, rho,
replicate, method, q, fdp, power, n_selected
```

Each simulation output directory contains scenario definitions, per-replicate
checkpoints, `replicate_metrics.csv`, `summary.csv`,
`method_diagnostics.csv`, and `run_manifest.rds`. A method failure is recorded
with its error message and represented by `NA`—never by a misleading zero FDP,
power, or selection count. Generated result directories are intentionally
ignored by Git.

## Runtime note

The full studies are computationally expensive, especially J=1000 because each
replicate generates feature-wise cross-fitted controls and evaluates 50 splits.
The code in this revision was validated with focused syntax, generator,
method-level, and plotting checks; the complete 200-replicate simulation grid
was not rerun.
