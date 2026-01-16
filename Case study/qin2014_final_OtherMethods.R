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
ncol(otu_filtered)
ncol(otu_t)

# ==============================================================================
# 3. Merge with phenotype (age) & basic setup
# ==============================================================================

names(metadata) 

# Check alignment
stopifnot(all(metadata$subject_id %in% rownames(otu_t)))

mydata <- merge(
  x     = metadata,
  y     = otu_filtered,
  by.x  = "subject_id",
  by.y  = "row.names"
)

# Drop SampleID after merge
mydata[, 1]
mydata <- mydata[, -1]
names(mydata)[1:35]

y <- mydata$albumine#log(mydata$Total_bilirubin)
hist(mydata$albumine)
names(mydata)
X=mydata[,-c(1:3,5,7,9:21,23,25:28)]
p=ncol(X)
names_x=names(X)
names(X)=paste0('X',1:p)
mydata_full=as.data.frame(cbind(y,X))
names(mydata_full)

# ==============================================================================
# 4. Parameters and  Setup
# ==============================================================================
amountTrain <- 0.5
amountTest  <- 1 - amountTrain
num_split   <- 50
q=0.1
set.seed(20260801)
# Quick fix
X$X1 <- as.numeric(as.factor(X$X1))-1
X$X3 <- as.numeric(as.factor(X$X3))-1
X$X5 <- as.numeric(as.factor(X$X5))-1

# Convert remaining to numeric (columns X3, X5, X7-X210)
X_matrix <- as.matrix(X)

# ==============================================================================
# 5. Dai data splitting
# ==============================================================================
Dai=function(X, y, num_split=5, q){
  n <- dim(X)[1]; p <- dim(X)[2]
  X_matrix <- as.matrix(X)
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  fdp <- rep(0, num_split)
  power <- rep(0, num_split)
  num_select <- rep(0, num_split)
  
  for(iter in 1:num_split){
    ### randomly split the data
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)
    
    ### get the penalty lambda for Lasso
    cvfit <- cv.glmnet(X_matrix[sample_index1, ], y[sample_index1], type.measure = "mse", nfolds = 10)
    lambda <- cvfit$lambda.min
    
    ### run Lasso on the first half of the data
    beta1 <- as.vector(glmnet(X_matrix[sample_index1, ], y[sample_index1], family = "gaussian", alpha = 1, lambda = lambda)$beta)
    nonzero_index <- which(beta1 != 0)
    if(length(nonzero_index)!=0){
      ### run OLS on the second half of the data, restricted on the selected features
      beta2 <- rep(0, p)
      beta2[nonzero_index] <- as.vector(lm(y[sample_index2] ~ as.matrix(X[sample_index2, nonzero_index]) - 1)$coeff)
      
      ### calculate the mirror statistics
      M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
      # M <- abs(beta1 + beta2) - abs(beta1 - beta2)
      selected_index <- analys(M, abs(M), q)
      
      ### number of selected variables
      if(length(selected_index)!=0){
        num_select[iter] <- length(selected_index)
        inclusion_rate[iter, selected_index] <- 1/num_select[iter]
        
        ### calculate fdp and power
        result <- fdp_power(selected_index, signal_index)
        fdp[iter] <- result$fdp
        power[iter] <- result$power
      }
    }
  }
  
  ### single data-splitting (DS) result
  DS_fdp <- fdp[1]
  DS_power <- power[1]
  
  ### multiple data-splitting (MDS) result
  inclusion_rate <- apply(inclusion_rate, 2, mean)
  
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
    return(selected_index)
}
Dai(X=X,y=y,num_split=num_split,q=0.1) 
# ==============================================================================
 # 5. BH
# ==============================================================================  
set.seed(20260801)
p <- dim(X)[2]
multi_fit <- multi.split(X_matrix, y, B = num_split)
pvalues = multi_fit$pval.corr
sorted_pvalues = sort(pvalues, decreasing = F, index.return = T)
cutoff <- which(sorted_pvalues$x <= (1:p)*q)

if(length(cutoff)){
  cutoff <- max(cutoff)
  selected_index <- sorted_pvalues$ix[1:cutoff]
}else{
  selected_index <- NA
}
names_x[selected_index]


# ==============================================================================
# 5. Knock-off
# ==============================================================================  

q = 0.1
nrounds = 100
param =list(
  objective = "reg:squarederror",
  eta       = 0.03,
  max_depth = 9,
  subsample = 0.8,
  colsample_bytree = 1,
  lambda    = 2,
  alpha     = 1
)
seed = 1292025

set.seed(seed)

## 1. build knockoff copies (second-order Gaussian)
ko <- create.second_order(X_matrix)             # returns list with X_k, Sigma, etc.
Xk <- ko

## 2. fit XGBoost on [X, X_k]
X_aug <- cbind(X_matrix, Xk)
colnames(X_aug) <- c(paste0("X",  seq_len(ncol(X))),
                     paste0("Xk", seq_len(ncol(X))))
dtrain <- xgb.DMatrix(data = X_aug, label = y)

bst <- xgb.train(params   = param,
                 data     = dtrain,
                 nrounds  = nrounds,
                 verbose  = 0)

## 3. feature-gain importance for originals & knockoffs
imp <- xgb.importance(model = bst, feature_names = colnames(X_aug))
gain <- numeric(ncol(X_aug)); names(gain) <- colnames(X_aug)
gain[imp$Feature] <- imp$Gain            # features not used keep gain = 0

## 4. antisymmetric statistic
p <- ncol(X)
W <- gain[1:p] - gain[(p + 1):(2 * p)]
names(W) <- paste0("X", seq_len(p))

## 5. apply knockoff threshold
T <- knockoff.threshold(W, fdr = q, offset = 1)   # BC+1 threshold
selected_index <- which(W >= T)
selected_index

length(selected_index)

names_x[selected_index]

 


