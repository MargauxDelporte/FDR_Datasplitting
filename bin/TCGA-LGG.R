setwd('C:/Users/mde4023/Downloads/TCGA-LGG')
BiocManager::install(c("TCGAbiolinks", "dplyr", "survival")) # survival is needed for outcome
library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
# --- 1. Define Queries for Features (X) ---
query_rnaseq <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

# You can add a second query here for Methylation if desired:
# query_meth <- GDCquery(...)

# --- 2. Download Data ---
GDCdownload(query_rnaseq) # Downloads the raw data

# --- 3. Prepare/Parse the Data ---
# Preprocessing step to harmonize and clean the data
rna_data <- GDCprepare(query = query_rnaseq, summarizedExperiment = TRUE)

# Extract the count matrix and perform variance filtering to reduce p
X_data_rnaseq <- as.data.frame(assay(rna_data))

# Filter X: Keep only the top 5000 most variable genes (adjust p here)
var_genes <- apply(X_data_rnaseq, 1, var)
X_data_rnaseq <- X_data_rnaseq[order(var_genes, decreasing=TRUE)[1:5000],]
X_data <- t(X_data_rnaseq) # Transpose to Samples (n) x Features (p)
# Extract Clinical Data
clinical <- GDCquery_clinic(project = "TCGA-LGG", type = "clinical")
names(clinical)
ANEUPLOIDY_SCORE
grep("ploid", colnames(colData(rna_data)), value = TRUE)
