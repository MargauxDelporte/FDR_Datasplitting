library(readr)
library(doParallel)
library(doRNG)      # Reproducible parallel random numbers
library(microbiome)
library(earth)
library(dplyr)
library(earth)        # MARS
library(foreach)
library(doParallel)
library(doRNG)  

source(paste0('C:/Users/mde4023/Downloads/FDR_Datasplitting','/Functions/HelperFunctions.R'))
# ==============================================================================
# ISSUEE: new levels (unseen values) bc only 48 samples in  set
# ==============================================================================

# Function to calculate Permuted R2 for a specific feature j
permR2Mars <- function(data, Y, j, model) {
  Xperm <- data
  # Permute column j
  Xperm[, j] <- sample(data[, j], replace = FALSE)
  
  # Predict using permuted data
  pred_perm <- predict(model, newdata = as.data.frame(Xperm))
  
  # Calculate R2
  rsq_perm <- 1 - sum((Y - pred_perm)^2) / sum((Y - mean(Y))^2)
  return(rsq_perm)
}

# Function to Select Features based on Mirror Statistics (FDR control)
# Recreates logic likely found in your original "SelectFeatures" function
SelectFeatures <- function(mm, ww, q){
  ### mm: mirror statistics mm=M
  ### ww: absolute value of mirror statistics ww=abs(M)
  ### q:  FDR control level q=0.1
  cutoff_set <- max(ww)
  for(t in ww){
    ps <- length(mm[mm > t])
    ng <- length(na.omit(mm[mm < -t]))
    rto <- (ng+1)/max(ps, 1)
    if(rto <= q){
      cutoff_set <- c(cutoff_set, t)
    }
  }
  cutoff <- min(cutoff_set)
  selected_index <- which(mm > cutoff)
  
  return(selected_index)
}

# ==============================================================================
# 3. Data Preparation
# ==============================================================================
setwd("C:/Users/mde4023/Downloads/MLRepo-master/MLRepo-master/datasets/qin2014")
task <- read_delim("mapping-orig.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)
otu    <- read_delim("otutable.txt", delim = "\t", 
                                 escape_double = FALSE, trim_ws = TRUE)
taxa   <- read_delim("taxatable.txt", delim = "\t", 
                     escape_double = FALSE, trim_ws = TRUE)


str(task)
dim(otu)
dim(taxa)
#But for modeling, it’s usually better to have samples as rows.
names=otu[,1]
otu_t <- t(otu[,-1]) 
# Prevalence of each OTU
prev <- colSums(otu_t > 0) / nrow(otu_t)
# Step 2: Compute relative abundances
otu_rel <- otu_t / rowSums(otu_t)

min_prev   <- 0.10     # present in at least 10% samples
min_relab  <- 1e-4     # relative abundance cutoff
min_relab_prev <- 0.05 # and observed in ≥5% samples

keep_prev <- prev >= min_prev
keep_abund <- colSums(otu_rel > min_relab) / nrow(otu_rel) >= min_relab_prev

keep <- keep_prev & keep_abund

otu_filtered <- otu_t[, keep]
row.names(otu_filtered)

names(task)[1]='SampleID'
#name task dataset
all(task$SampleID %in%rownames(otu_filtered))


#merge task and otu datasets
mydata <- merge(task, otu_filtered, by.x="SampleID", by.y="row.names")
mydata=mydata[,-1]
#View(data)
names(mydata)
hist(mydata$Total_bilirubin, main = "Total_bilirubin Distribution", xlab = "Total_bilirubin")
hist(log(mydata$Total_bilirubin), main = "Log Total_bilirubin Distribution", xlab = "Total_bilirubin")

# ==============================================================================
# 4. Parameters and Parallel Setup
# ==============================================================================
amountTrain <- 0.5
amountTest  <- 1 - amountTrain
num_split   <- 10   # Number of data splits
q_level     <- 0.10 # Target FDR level

# Setup Parallel Backend
n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1)
cl <- parallel::makeCluster(n_cores)
registerDoParallel(cl)
registerDoRNG(123) # Set seed for reproducibility
names(data)
y=log(mydata$Total_bilirubin)
X=mydata[,20:ncol(mydata)]
# names(X)=paste0('X',1:p)

data=cbind(y,X)
n= nrow(mydata)
p=ncol(mydata)-1
# ==============================================================================
# 5. Main Parallel Loop (Data Splitting)
# ==============================================================================
res_mat <- foreach(iter = 1:num_split,
                   .combine = "rbind",
                   .packages = c("earth"),
                   .export   = c("permR2Mars","SelectFeatures","CalculateFDP_Power")) %dorng% {
                     
                     # --- indices ---
                     train_index <- sample.int(n, size = floor(amountTrain * n), replace = FALSE)
                     remaining_index <- setdiff(seq_len(n), train_index)

                     # split the remaining half evenly
                     size_half <- floor((amountTest/2) * n)
                     sample_index1 <- sample(remaining_index, size = size_half, replace = FALSE)
                     sample_index2 <- setdiff(remaining_index, sample_index1)
                     dataTrain <- data[train_index, , drop = FALSE]
                     
                     mars_poly= earth(
                       y ~ .,
                        data    = dataTrain
                     )
                     lm <- mars_poly
                     lm
                     # --- R^2 on train (not stored) / test halves (stored) ---?earth
                     pred1 <- predict(lm, newdata = X[sample_index1,])
                     pred2 <- predict(lm, newdata = X[sample_index2,])
                     y1 <- y[sample_index1]; y2 <- y[sample_index2]
                     
                     R2orig1 <- 1 - sum((y1 - pred1)^2) / sum((y1 - mean(y1))^2)
                     R2orig2 <- 1 - sum((y2 - pred2)^2) / sum((y2 - mean(y2))^2)
                     R2orig1;R2orig2
                     # --- permutation-based drops ---
                     Rnew1 <- sapply(seq_len(p), function(j)
                       permR2Mars(as.data.frame(X[sample_index1, , drop = FALSE]), Y = y1, j = j, model = lm))
                     Rnew2 <- sapply(seq_len(p), function(j)
                       permR2Mars(as.data.frame(X[sample_index2, , drop = FALSE]), Y = y2, j = j, model = lm))
                     
                     beta1 <- R2orig1 - Rnew1
                     beta2 <- R2orig2 - Rnew2
                     mirror <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
                     #hist(mirror)
                     selected_index <- SelectFeatures(mirror, abs(mirror),q=0.1)
                     num_sel <- length(selected_index)
                     num_sel
                     inc_row <- numeric(p)
                     fdp_val <- 0
                     pow_val <- 0
                     
                     if (num_sel > 0) {
                       inc_row[selected_index] <- 1 / num_sel
                     }
                     c(num_sel, R2orig1, R2orig2, inc_row)
                   }

parallel::stopCluster(cl)

# ---- unpack ----
num_select     <- res_mat[, 1]
R2orig1_vec    <- res_mat[, 2]
R2orig2_vec    <- res_mat[, 3]
inclusion_rate_mat <- res_mat[, -(1:3), drop = FALSE]  # num_split x p

# MDS inclusion rates (avg across splits)
inclusion_rate <- apply(inclusion_rate_mat, 2, mean)

# Rank features by empirical inclusion rate
feature_rank <- order(inclusion_rate)
feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))

if (length(feature_rank) != 0) {
  null_feature <- numeric()
  for (feature_index in seq_along(feature_rank)) {
    if (sum(inclusion_rate[feature_rank[1:feature_index]]) > q) break
    null_feature <- c(null_feature, feature_rank[feature_index])
  }
  selected_index <- setdiff(feature_rank, null_feature)
}
selected_index
mean(R2orig1_vec)
mean(R2orig2_vec)
