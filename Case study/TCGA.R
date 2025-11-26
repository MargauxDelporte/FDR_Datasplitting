## ---- Install packages (run once) ----
# install.packages(c("BiocManager", "glmnet", "dplyr"))
# BiocManager::install("TCGAbiolinks")
# BiocManager::install("immunedeconv")
## 1. Set up a personal library (change username if needed)
#dir.create("C:/Users/mde4023/Rlibs", showWarnings = FALSE)
.libPaths("C:/Users/mde4023/Rlibs")
.libPaths()   # check that this path is now first

setwd('C:/Users/mde4023/OneDrive - Weill Cornell Medicine/Documents/TCGA_downloads')
## 2. Install remotes (into your user library, not Program Files)
##install.packages("remotes")

## 3. Install immunedeconv from GitHub (this pulls CRAN + Bioconductor deps)
##remotes::install_github("omnideconv/immunedeconv")
library(biomaRt)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(immunedeconv)
library(glmnet)
library(dplyr)
library(TCGAbiolinks)
project <- "TCGA-BRCA"
my_dir='C:/Users/mde4023/OneDrive - Weill Cornell Medicine/Documents/TCGA_downloads'

query <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
Sys.setenv("GDC_CLIENT_DIR" = "C:/Users/mde4023/OneDrive - Weill Cornell Medicine/Documents/TCGA_downloads/gdc-client")
GDCdownload(
  query,
  directory = my_dir,
  method    = "client",
  files.per.chunk = 5
)

brca_expr <- GDCprepare(query, directory = my_dir)

## brca_expr is a SummarizedExperiment:
# 1. Extract raw counts matrix
expr_mat <- assay(brca_expr)   # genes x samples

# 2. Clean Ensembl IDs: remove version (e.g. ENSG00000141510.14 → ENSG00000141510)
ensembl_ids <- rownames(expr_mat)
ensembl_clean <- sub("\\..*$", "", ensembl_ids)

# 3. Query Ensembl (GRCh38)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

annot <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique(ensembl_clean),
  mart = mart
)

# 4. Merge annotation with expression matrix
annot_df <- data.frame(
  ensembl_gene_id = ensembl_clean,
  row_id = seq_along(ensembl_clean),
  stringsAsFactors = FALSE
)

annot_merged <- left_join(annot_df, annot, by = "ensembl_gene_id")

# 5. Replace empty symbols with NA
annot_merged$hgnc_symbol[annot_merged$hgnc_symbol == ""] <- NA

# 6. Keep only genes with valid gene symbols
keep_idx <- which(!is.na(annot_merged$hgnc_symbol))

expr_kept <- expr_mat[keep_idx, ]
symbols_kept <- annot_merged$hgnc_symbol[keep_idx]

# 7. Aggregate duplicates by mean
expr_df <- as.data.frame(expr_kept)
expr_df$symbol <- symbols_kept

expr_sym_agg <- expr_df %>%
  group_by(symbol) %>%
  summarize(across(where(is.numeric), mean), .groups = "drop")

# 8. Final expression matrix: gene symbols × samples
expr_final <- as.matrix(expr_sym_agg[, -1])
rownames(expr_final) <- expr_sym_agg$symbol
# Now rows = gene symbols, cols = samples
dim(expr_final)

## ---- 3. Compute immune cell infiltration scores with immunedeconv ----
# Set default deconvolution method
set_cibersort_binary("CIBERSORT.R")  # path to CIBERSORT script if you use it
# Alternatively use xCell, quanTIseq, etc.
deconvolute_func <- deconvolute_xcell  # for example

# immunedeconv expects a data.frame with genes in first column, samples as others
expr_for_deconv <- data.frame(gene = rownames(expr_final),
                              expr_final,
                              check.names = FALSE)

immune_res <- deconvolute_func(expr_for_deconv)
?deconvolute_func

# immune_res: cell types x samples; first column is 'cell_type'
immune_mat <- immune_res[, -1]
rownames(immune_mat) <- immune_res$cell_type

# Choose an outcome, e.g. "CD8+ T-cells" (name depends on method)
rownames(immune_mat)
# Suppose the row is named "CD8+ T-cells" (adjust string if needed)
y <- as.numeric(immune_mat["CD8+ T-cells", ])

# Match samples / transpose predictors
X <- t(expr_final)  # samples x genes
X <- X[rownames(immune_mat), , drop = FALSE]  # align samples if needed

## ---- 4. Train / test split ----
set.seed(123)
n <- nrow(X)
train_idx <- sample(seq_len(n), size = floor(0.7 * n))

X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]
y_train <- y[train_idx]
y_test  <- y[-train_idx]

## ---- 5. Fit elastic net (glmnet) ----
# Standardize predictors inside glmnet
alpha_val <- 0.5  # elastic net
cvfit <- cv.glmnet(
  x = X_train,
  y = y_train,
  alpha = alpha_val,
  family = "gaussian",
  nfolds = 5,
  parallel = FALSE
)

best_lambda <- cvfit$lambda.min
best_lambda

# Fit final model
fit <- glmnet(
  x = X_train,
  y = y_train,
  alpha = alpha_val,
  lambda = best_lambda,
  family = "gaussian"
)

## ---- 6. Evaluate R² on test set ----
pred_test <- predict(fit, newx = X_test, s = best_lambda)[, 1]

rss  <- sum((y_test - pred_test)^2)
tss  <- sum((y_test - mean(y_test))^2)
R2   <- 1 - rss / tss
R2                
