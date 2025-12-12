library(readr)
library(doParallel)
library(doRNG)      # Reproducible parallel random numbers
library(microbiome)
library(earth)
library(dplyr)
library(earth)        # MARS
library(foreach)
library(doParallel)
library(ranger)  
library(randomForest)
source(paste0('C:/Users/mde4023/Downloads/FDR_Datasplitting','/Functions/HelperFunctions.R'))
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

names(task) 

# Check alignment
stopifnot(all(task$SampleID %in% rownames(otu_t)))

mydata <- merge(
  x     = task,
  y     = otu_filtered,
  by.x  = "SampleID",
  by.y  = "row.names"
)

# Drop SampleID after merge
mydata[, 1]
mydata <- mydata[, -1]

mydata=mydata[complete.cases(mydata),]
n <- nrow(mydata)
p <- ncol(mydata) - 1   # number of OTUs

y <- log(mydata$Total_bilirubin)
names(mydata)
#X=mydata[,c(1:19)]
X=mydata[,-c(1,4,8,9,12,13,16,17,18,19)]


# ==============================================================================
# 4. Parameters and Parallel Setup
# ==============================================================================
amountTrain <- 0.5
amountTest  <- 1 - amountTrain
num_split   <- 50#5010   # Number of data splits
q=0.1
# Setup Parallel Backend
#X=mydata[,-c(1,4,8,9,12,13,16,17,18,19)]
p <- ncol(X)   # number of OTUs
n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1)
cl <- parallel::makeCluster(n_cores)
registerDoParallel(cl)
registerDoRNG(11272025) # Set seed for reproducibility
# ==============================================================================
# 5. Main Parallel Loop (Data Splitting)
# ==============================================================================
res_mat <- foreach(iter = 1:num_split,
                   .combine = "rbind",
                   .packages = c("randomForest")) %dorng% {
                     # --- indices ---
                     train_index <- sample.int(n, size = floor(amountTrain * n), replace = FALSE)
                     remaining_index <- setdiff(seq_len(n), train_index)
                     data=cbind(y,X)
                     # split the remaining half evenly ?earth
                     size_half <- floor((amountTest/2) * n)
                     sample_index1 <- sample(remaining_index, size = size_half, replace = FALSE)
                     sample_index2 <- setdiff(remaining_index, sample_index1)
                     dataTrain <- data[train_index, , drop = FALSE]
                     mars_poly=  randomForest(
                       y ~ ., 
                       ntry=100,
                       ntree=104,
                       nodesize=1,
                       data    = data_full
                     )
                     
                     
                     lm <- mars_poly   # keep the rest of your code unchanged
                     
                     lm
                     
                     ## --- R² on the two halves, same as before ---
                     
                     pred1 <- predict(lm, newdata = as.data.frame(X[sample_index1, ]))
                     pred2 <- predict(lm, newdata = as.data.frame(X[sample_index2, ]))
                     
                     y1 <- y[sample_index1]
                     y2 <- y[sample_index2]
                     
                     R2orig1 <- 1 - sum((y1 - pred1)^2) / sum((y1 - mean(y1))^2)
                     R2orig2 <- 1 - sum((y2 - pred2)^2) / sum((y2 - mean(y2))^2)
                     
                     R2orig1
                     R2orig2
                     # --- permutation-based drops ---
                     Rnew1 <- sapply(seq_len(p), function(j)
                       permR2(as.data.frame(X[sample_index1, , drop = FALSE]), Y = y1, j = j, model = lm))
                     Rnew2 <- sapply(seq_len(p), function(j)
                       permR2(as.data.frame(X[sample_index2, , drop = FALSE]), Y = y2, j = j, model = lm))
                     
                     beta1 <- R2orig1 - Rnew1
                     beta2 <- R2orig2 - Rnew2
                     mirror <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
                     #hist(mirror)
                     #summary(mirror)
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
selected_index
mean(R2orig1_vec)
mean(R2orig2_vec)
