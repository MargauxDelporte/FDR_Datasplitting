### High dimension linear model
rm(list = ls())

mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
setwd(mywd)

source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Functions/TriangleBoosterTrainMS.R'))
source(paste0(mywd,'/Functions/ApplyGBMKnockoff.R'))
source(paste0(mywd,'/Functions/MarsParallelHD_new.R'))

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




# ==============================================================================
# 1. File Paths & Data Import
# ==============================================================================

base_dir <- "C:/Users/mde4023/Downloads/FDR_Datasplitting/Case study/qin2014"

task <- read_delim(file.path(base_dir, "mapping-orig.txt"), 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)
names(task)[1]='SampleID'
otu <- read_delim(
  file         = file.path(base_dir, "otutable.txt"),
  delim        = "\t",
  escape_double = FALSE,
  trim_ws       = TRUE
)

taxa <- read_delim(
  file         = file.path(base_dir, "taxatable.txt"),
  delim        = "\t",
  escape_double = FALSE,
  trim_ws       = TRUE
)


# ==============================================================================
# 2. OTU Pre-processing (prevalence + relative abundance filter)
# ==============================================================================

# Samples as rows, OTUs as columns
otu_t <- t(as.matrix(otu[, -1]))   # drop OTU ID column, transpose
namess=as.vector(otu[,1])
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

namess_filtered=unlist(namess)[which(keep)]
# ==============================================================================
# 3. Merge with phenotype (age) & basic setup
# ==============================================================================
mydata <- merge(
  x     = task,
  y     = otu_filtered,
  by.x  = "SampleID",
  by.y  = "row.names"
)

# Drop SampleID after merge
mydata <- mydata[, -1]

mydata=mydata[complete.cases(mydata),]
n <- nrow(mydata)
p <- ncol(mydata) - 1   # number of OTUs

y <- mydata$Serum_albumin#log(mydata$Total_bilirubin)
X=mydata[,-c(1,4,8,9,11,12,13,16,17,18,19)]

namess_filtered_combined=c(names(task)[-1],namess_filtered)
namess_filtered_combined=namess_filtered_combined[-c(1,4,8,9,11,12,13,16,17,18,19)]
mydata_full=as.data.frame(cbind(y,X))
p=1042;n=130
data_full=mydata_full
X=mydata_full[,-1]
y=mydata_full[,1]
#table(X[,1])
for(i in 1:8){
  print(names(X)[i])
  print(table(X[,i]))
}
table(X$Ascites)
X$Sex=ifelse(X$Sex=='female',1,0)
X$Cirrhotic=ifelse(X$Cirrhotic=='Cirrhosis',1,0)
X$HBV_related=ifelse(X$HBV_related=='Y',1,0)
X$Alcohol_related=ifelse(X$Alcohol_related=='Y',1,0)
X$Ascites=ifelse(X$Ascites=='Mild'|X$Ascites=='Sever',1,0)
X$HE=ifelse(X$HE=='Grade 1',1,0)
# ==============================================================================
# 4. Parameters and Parallel Setup
# ==============================================================================
q=0.1
num_split=5
set.seed(1292025)
# ==============================================================================
# 5. Dai data splitting
# ==============================================================================
Dai=function(X, y, num_split=5, q){
  n <- dim(X)[1]; p <- dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  fdp <- rep(0, num_split)
  power <- rep(0, num_split)
  num_select <- rep(0, num_split)
  
  for(iter in 1:num_split){
    ### randomly split the data
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)
    
    ### get the penalty lambda for Lasso
    cvfit <- cv.glmnet(X[sample_index1, ], y[sample_index1], type.measure = "mse", nfolds = 10)
    lambda <- cvfit$lambda.min
    
    ### run Lasso on the first half of the data
    beta1 <- as.vector(glmnet(X[sample_index1, ], y[sample_index1], family = "gaussian", alpha = 1, lambda = lambda)$beta)
    nonzero_index <- which(beta1 != 0)
    if(length(nonzero_index)!=0){
      ### run OLS on the second half of the data, restricted on the selected features
      beta2 <- rep(0, p)
      beta2[nonzero_index] <- as.vector(lm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1)$coeff)
      
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
X_df  <- as.data.frame(X)           # list → data.frame
X_mat <- as.matrix(X_df)            # data.frame → numeric matrix
Dai(X=X_mat,y=y,num_split=5,q=0.1) 

# ==============================================================================
 # 5. BH
# ==============================================================================  
set.seed(1292025)
p <- dim(X)[2]
multi_fit = multi.split(X_mat, y, B = num_split)
pvalues = multi_fit$pval.corr
sorted_pvalues = sort(pvalues, decreasing = F, index.return = T)
cutoff <- which(sorted_pvalues$x <= (1:p)*q)

if(length(cutoff)){
  cutoff <- max(cutoff)
  selected_index <- sorted_pvalues$ix[1:cutoff]
}else{
  selected_index <- NA
}
names(X)[selected_index]


# ==============================================================================
# 5. Knock-off
# ==============================================================================  

q = 0.1
nrounds = 500
param =list(
  objective = "reg:squarederror",
  eta       = 0.05,
  max_depth = 3,
  subsample = 0.6,
  colsample_bytree = 0.8,
  lambda    = 1,
  alpha     = 0
)
seed = 1292025

set.seed(seed)

## 1. build knockoff copies (second-order Gaussian)
ko <- create.second_order(X)             # returns list with X_k, Sigma, etc.
Xk <- ko

## 2. fit XGBoost on [X, X_k]
X_aug <- cbind(X, Xk)
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
#X3   X7  X21  X30  X57 X186 X442 X516 X574 X584 X585 X630 X916 
#3    7   21   30   57  186  442  516  574  584  585  630  916 


#[1] if num_split=5
#[1] 574 584 585 592 611 615 630 683 882 896
length(selected_index)

namess_filtered_combined[selected_index]


GBM_knockoff(X=X_mat,y=y,q=0.1) 


