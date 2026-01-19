rm(list = ls())

mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
setwd(mywd)
source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Functions/ApplyGBMKnockoff.R'))


source(paste0(mywd,'/Functions Dai/knockoff.R'))
source(paste0(mywd,'/Functions Dai/analysis.R'))
source(paste0(mywd,'/Functions Dai/MBHq.R'))
source(paste0(mywd,'/Functions Dai/DS.R'))
source(paste0(mywd,'/Functions Dai/fdp_power.R'))

#devtools::install_github("Jeremy690/DSfdr/DSfdr",force = TRUE)
library(xgboost)
library(gbm)
library(ranger)
library(MASS)

library(glmnet)
library(knockoff)
library(mvtnorm)
library(hdi)

library(readr)
library(doParallel)
library(doRNG)      # Reproducible parallel random numbers
library(dplyr)
library(foreach)

# Load library
library(curatedMetagenomicData)
library(dplyr)

# ==============================================================================
# 1. File Paths & Data Import
# ==============================================================================
# Get the data - returns a TreeSummarizedExperiment object
qin_data <- curatedMetagenomicData(
  "QinN_2014.relative_abundance",
  dryrun = FALSE,
  counts = TRUE
)


# Access the data
qin_tse <- qin_data[[1]]  # Extract TreeSummarizedExperiment

# Get relative abundance matrix
abundance <- assay(qin_tse)

# Get sample metadata
metadata <- colData(qin_tse)
table(metadata$disease)
metadata$ascites <- ifelse(grepl("ascites", metadata$disease), 1, 0)
metadata$hepatitis <- ifelse(grepl("hepatitis", metadata$disease), 1, 0)
metadata$schistosoma <- ifelse(grepl("schistosoma", metadata$disease), 1, 0)
metadata$cirrhosis <- ifelse(grepl("cirrhosis", metadata$disease), 1, 0)

# Get taxonomic information
taxonomy <- rowData(qin_tse)

# ==============================================================================
#Helper functions
# ==============================================================================

# Function to calculate Permuted R2 for a specific feature j
permR2 <- function(data, Y, j, model) {
  Xperm <- data
  # Permute column j
  Xperm[, j] <- sample(data[, j], replace = FALSE)
  
  # Predict using permuted data
  pred_perm <- predict(model, newdata = as.data.frame(Xperm))
  
  # Calculate R2
  rsq_perm <- 1 - sum((Y - pred_perm)^2) / sum((Y - mean(Y))^2)
  return(rsq_perm)
}


# ==============================================================================
# 2. OTU Pre-processing (prevalence + relative abundance filter)
# ==============================================================================

# Samples as rows, OTUs as columns
otu_t <- t(abundance)   # drop OTU ID column, transpose
namess=row.names(abundance)
# Prevalence per OTU
prev <- colSums(otu_t > 0) / nrow(otu_t)

# Relative abundances
otu_rel <- otu_t / rowSums(otu_t)

# Filtering thresholds
min_prev        <- 0.10   # present in ≥10% of samples
min_relab       <- 1e-4   # minimum relative abundance
min_relab_prev  <- 0.05   # abundance > min_relab in ≥5% of samples

keep_prev  <- prev >= min_prev
keep_abund <- colSums(otu_rel > min_relab) / nrow(otu_rel) >= min_relab_prev
keep       <- keep_prev & keep_abund

otu_filtered <- otu_t[, keep, drop = FALSE]

# ==============================================================================
# 3. Merge with phenotype (age) & basic setup
# ==============================================================================

# Check alignment
stopifnot(all(metadata$subject_id %in% rownames(otu_t)))

mydata <- merge(
  x     = metadata,
  y     = otu_filtered,
  by.x  = "subject_id",
  by.y  = "row.names"
)

# Drop SampleID and irrelevant predictors after merge
mydata[, 1]
mydata <- mydata[, -1]
y <- mydata$albumine
X=mydata[,-c(1:5,7,9:21,23,25:28)]
p=ncol(X)

#rename the columns
names_x=names(X)
names(X)=paste0('X',1:p)
mydata_full=as.data.frame(cbind(y,X))

# ==============================================================================
# 4. Parameters and  Setup
# ==============================================================================
amountTrain <- 0.5
amountTest  <- 1 - amountTrain
num_split   <- 50
q=0.1
set.seed(20260801)
# Quick fix
names_x[1:8]
X$X2 <- as.numeric(as.factor(X$X2))-1
X$X4 <- as.numeric(as.factor(X$X4))-1
n <- dim(X)[1]; p <- dim(X)[2]
# Convert remaining to numeric (columns X3, X5, X7-X210)
X_matrix <- as.matrix(X)

# ==============================================================================
#Helper functions
# ==============================================================================

# Function to calculate Permuted R2 for a specific feature j
permR2 <- function(data, Y, j, model) {
  Xperm <- data
  # Permute column j
  Xperm[, j] <- sample(data[, j], replace = FALSE)
  
  # Predict using permuted data
  pred_perm <- predict(model, newdata = as.data.frame(Xperm))
  
  # Calculate R2
  rsq_perm <- 1 - sum((Y - pred_perm)^2) / sum((Y - mean(Y))^2)
  return(rsq_perm)
}



detectCores()
cl <- parallel::makeCluster(10)
registerDoParallel(cl)

registerDoRNG(20261901) 
set.seed(20261901)
# ==============================================================================
# 5. Main Parallel Loop (Data Splitting)
# ==============================================================================
param = list(
  objective = "reg:squarederror",
  eta       = 0.010,
  max_depth = 7,
  subsample = 0.6,
  colsample_bytree = 0.6,
  lambda    = 0,
  alpha     = 0.5
)
nrounds = 100

res_mat <- foreach(iter = 1:num_split,
                   .combine = "rbind",
                   .packages = c("xgboost")) %dorng% {
                     set.seed(iter)
                     permR2 <- function(data, Y, j, model) {
                       Xperm <- data
                       # Permute column j
                       Xperm[, j] <- sample(data[, j], replace = FALSE)
                       
                       # Predict using permuted data
                       pred_perm <- predict(model, newdata = as.data.frame(Xperm))
                       
                       # Calculate R2
                       rsq_perm <- 1 - sum((Y - pred_perm)^2) / sum((Y - mean(Y))^2)
                       return(rsq_perm)
                     }
                     source(paste0('C:/Users/mde4023/Downloads/FDR_Datasplitting','/Functions/HelperFunctions.R'))
                     data_full=mydata_full
                     X=X_matrix
                                # --- indices ---
                     train_index     <- sample.int(n, size = floor(amountTrain * n), replace = FALSE)
                     remaining_index <- setdiff(seq_len(n), train_index)
                     
                     # split the remaining part in two halves
                     size_half     <- floor((amountTest / 2) * n)
                     sample_index1 <- sample(remaining_index, size = size_half, replace = FALSE)
                     sample_index2 <- setdiff(remaining_index, sample_index1)
                     
                     # --- fit RF using parameter vector pm ---350 260 16
                     dtrain <- xgb.DMatrix(data = X[train_index,], label = y[train_index])
                     bst <- xgb.train(params   = param,
                                      data     = dtrain,
                                      nrounds  = nrounds,
                                      verbose  = 0)
                     
                     # --- R² on the two halves ---
                     pred1 <- predict(bst, newdata = as.data.frame(X[sample_index1, , drop = FALSE]))
                     pred2 <- predict(bst, newdata = as.data.frame(X[sample_index2, , drop = FALSE]))
                     
                     y1 <- y[sample_index1]
                     y2 <- y[sample_index2]
                     
                     R2orig1 <- 1 - sum((y1 - pred1)^2) / sum((y1 - mean(y1))^2)
                     R2orig2 <- 1 - sum((y2 - pred2)^2) / sum((y2 - mean(y2))^2)
                     
                     # --- permutation-based drops ---
                     Rnew1 <- sapply(seq_len(p), function(j)
                       permR2(as.data.frame(X[sample_index1, , drop = FALSE]), Y = y1, j = j, model = bst))
                     Rnew2 <- sapply(seq_len(p), function(j)
                       permR2(as.data.frame(X[sample_index2, , drop = FALSE]), Y = y2, j = j, model = bst))
                     
                     beta1  <- R2orig1 - Rnew1
                     beta2  <- R2orig2 - Rnew2
                     mirror <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
                     selected_index <- SelectFeatures(mirror, abs(mirror), q = 0.10)
                     num_sel <- length(selected_index)
                     
                     inc_row <- numeric(p)
                     if (num_sel > 0) {
                       inc_row[selected_index] <- 1 / num_sel
                     }
                     
                     c(num_sel, R2orig1, R2orig2, inc_row)
                   }

# ---- unpack ----
num_select     <- res_mat[, 1]
R2orig1_vec    <- res_mat[, 2]
R2orig2_vec    <- res_mat[, 3]

inclusion_rate_mat <- res_mat[, -(1:3), drop = FALSE]  # num_split x p
mean(c(R2orig1_vec,R2orig2_vec))
# MDS inclusion rates (avg across splits)
### multiple data-splitting (MDS) result
inclusion_rate <- apply(inclusion_rate_mat, 2, mean)

### rank the features by the empirical inclusion rate
feature_rank <- order(inclusion_rate)
feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
if(length(feature_rank)!=0){
  null_feature <- numeric()
  
  ### backtracking 
  for(feature_index in 1:length(feature_rank)){
    if(sum(inclusion_rate[feature_rank[1:feature_index]]) > q){
      break
    }else{
      null_feature <- c(null_feature, feature_rank[feature_index])
    }
  }
  selected_index <- setdiff(feature_rank, null_feature)
}

names_x[selected_index]
