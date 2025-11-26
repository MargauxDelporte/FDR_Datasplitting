library(doParallel)
library(doRNG)      # Reproducible parallel random numbers
BiocManager::install("microbiome")
library(microbiome)
library(earth)
library(dplyr)
library(earth)        # MARS
library(foreach)
library(doParallel)
library(doRNG)  

# ==============================================================================
# 2. Helper Functions
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
data("dietswap")
ps <- dietswap

# Extract Outcome (Shannon Diversity)
outcome <- estimate_richness(ps)$Shannon

# Extract Features (Centered Log Ratio transformed OTUs)
# Note: Ensure otu_table is oriented correctly (Taxa as rows initially)
otu_clr <- as.data.frame(t(otu_table(ps))) 

# Basic renaming for clarity
y <- outcome
X <- otu_clr
n <- nrow(X)
p <- ncol(X)
colnames(X) <- paste0("X", 1:p)

# Full Data Frame
df_full <- data.frame(y = y, X)

# Sanity Check: Initial Fit
# Note: Removed 's="lambda.min"' as that is for glmnet, not earth
fit_init <- earth(x = as.matrix(X), y = y)
pred_init <- predict(fit_init, newdata = as.matrix(X))
cat("Initial Model R2:", cor(y, pred_init)^2, "\n")

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

# ==============================================================================
# 5. Main Parallel Loop (Data Splitting)
# ==============================================================================

res_mat <- foreach(iter = 1:num_split, 
                   .combine = "rbind", 
                   .packages = c("earth"),
                   .export = c("permR2Mars", "select_features_fdr")) %dopar% {
                     
                     # --- A. Data Splitting ---
                     train_idx <- sample.int(n, size = floor(amountTrain * n), replace = FALSE)
                     remaining_idx <- setdiff(seq_len(n), train_idx)
                     
                     # Split remainder into two inference halves
                     size_half <- floor((amountTest / 2) * n)
                     sample_idx1 <- sample(remaining_idx, size = size_half, replace = FALSE)
                     sample_idx2 <- setdiff(remaining_idx, sample_idx1)
                     
                     dataTrain <- df_full[train_idx, , drop = FALSE]
                     
                     # --- B. Model Fitting (on Train) ---
                     mars_model <- earth(y ~ ., data = dataTrain, pmethod = "seqrep", thresh = 0)
                     
                     # --- C. Prediction on Inference Halves ---
                     # Helper to get R2
                     calc_r2 <- function(model, x_idx, y_vec) {
                       preds <- predict(model, newdata = df_full[x_idx, -1, drop=FALSE]) # -1 to exclude y
                       1 - sum((y_vec - preds)^2) / sum((y_vec - mean(y_vec))^2)
                     }
                     
                     y1 <- y[sample_idx1]
                     y2 <- y[sample_idx2]
                     
                     R2orig1 <- calc_r2(mars_model, sample_idx1, y1)
                     R2orig2 <- calc_r2(mars_model, sample_idx2, y2)
                     
                     # --- D. Permutation Importance (Mirror Statistics) ---
                     # Calculate importance on Half 1
                     Rnew1 <- sapply(seq_len(p), function(j) {
                       permR2Mars(df_full[sample_idx1, -1, drop=FALSE], Y = y1, j = j, model = mars_model)
                     })
                     
                     # Calculate importance on Half 2
                     Rnew2 <- sapply(seq_len(p), function(j) {
                       permR2Mars(df_full[sample_idx2, -1, drop=FALSE], Y = y2, j = j, model = mars_model)
                     })
                     
                     # Calculate Mirror Statistics
                     # Importance = Reduction in R2 (Original - Permuted)
                     beta1 <- R2orig1 - Rnew1
                     beta2 <- R2orig2 - Rnew2
                     
                     # Mirror stat: sign is product of signs, magnitude is sum of absolute importance
                     mirror <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
                     
                     # --- E. Feature Selection ---
                     selected_index <- select_features_fdr(mirror, q = q_level)
                     num_sel <- length(selected_index)
                     
                     # Create inclusion vector (1/num_sel for selected features)
                     inc_row <- numeric(p)
                     if (num_sel > 0) {
                       inc_row[selected_index] <- 1 / num_sel
                     }
                     
                     # Return: [NumSelected, R2_1, R2_2, InclusionVector...]
                     c(num_sel, R2orig1, R2orig2, inc_row)
                   }

# Stop cluster
parallel::stopCluster(cl)

# ==============================================================================
# 6. Results Aggregation
# ==============================================================================
# Unpack results
num_select_vec <- res_mat[, 1]
R2orig1_vec    <- res_mat[, 2]
R2orig2_vec    <- res_mat[, 3]
inclusion_rate_mat <- res_mat[, -(1:3), drop = FALSE]

# Calculate Average Inclusion Rate (Stability)
inclusion_rate <- colMeans(inclusion_rate_mat)

# Rank features
feature_rank <- order(inclusion_rate, decreasing = FALSE) # Lowest to Highest for loop check
feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))

# Final Selection Aggregation (heuristic based on inclusion sum)
final_selected_index <- integer(0)

if (length(feature_rank) > 0) {
  # Reverse to process most important first
  feature_rank_desc <- rev(feature_rank)
  null_features <- numeric()
  
  # Note: This logic assumes 'feature_rank' was intended to be checked against 'q'
  # Interpreting your original loop: it accumulates features until sum > q? 
  # Usually, for aggregation, we just take features with inclusion_rate > threshold (e.g. 0.5)
  # But sticking to your specific logic:
  
  for (idx in seq_along(feature_rank)) {
    # Check cumulative inclusion of the lowest ranked features
    if (sum(inclusion_rate[feature_rank[1:idx]]) > q_level) break
    null_features <- c(null_features, feature_rank[idx])
  }
  
  final_selected_index <- setdiff(feature_rank, null_features)
}

# --- Outputs ---
cat("\nAnalysis Complete.\n")
cat("Average R2 (Half 1):", mean(R2orig1_vec), "\n")
cat("Average R2 (Half 2):", mean(R2orig2_vec), "\n")
cat("Indices of Selected Features:", final_selected_index, "\n")
cat("Names of Selected Features:", colnames(X)[final_selected_index], "\n")

