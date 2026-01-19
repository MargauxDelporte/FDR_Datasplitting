library(xgboost)
mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting/Case Study'
setwd(mywd)
X_matrix <- read.csv('X_matrix_Qin2014.csv')
y <- read.csv('y.csv')
y <- y[,1]  # Convert to vector if it's a data frame
source('HelperFunctions.R')
# Function to calculate R²
calculate_r2 <- function(actual, predicted) {
  1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
}

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


# Define hyperparameter grid
param_grid <- expand.grid(
  eta = c(0.01, 0.05, 0.1, 0.3),
  max_depth = c(3, 5, 7, 9),
  subsample = c(0.6, 0.8, 1.0),
  colsample_bytree = c(0.6, 0.8, 1.0),
  lambda = c(0, 0.5, 1, 2),
  alpha = c(0, 0.5, 1),
  nrounds = c(100, 200, 500)
)

# Initialize results storage
results <- data.frame(
  eta = numeric(),
  max_depth = integer(),
  subsample = numeric(),
  colsample_bytree = numeric(),
  lambda = numeric(),
  alpha = numeric(),
  nrounds = integer(),
  R2_test1 = numeric(),
  R2_test2 = numeric(),
  R2_avg = numeric(),
  selected= numeric()
)

n <- nrow(X_matrix)
p <- ncol(X_matrix)

# --- indices ---
set.seed(20260501)
amountTrain <- 0.5
amountTest  <- 1 - amountTrain

train_index     <- sample.int(n, size = floor(0.5 * n), replace = FALSE)
remaining_index <- setdiff(seq_len(n), train_index)

# split the remaining part in two halves
size_half     <- floor((amountTest / 2) * n)
sample_index1 <- sample(remaining_index, size = size_half, replace = FALSE)
sample_index2 <- setdiff(remaining_index, sample_index1)

# Prepare data for XGBoost (requires matrix format)
X_train <- as.matrix(X_matrix[train_index, ])
y_train <- y[train_index]

# Grid search i=1
best_r2 <- -Inf
best_params <- NULL
most_selected=0
cat("Starting hyperparameter tuning...\n")
cat("Total combinations to test:", nrow(param_grid), "\n\n")

for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  
  # Train model with current parameters
  mynlm <- xgboost(
    x = X_train,
    y = y_train,
    max_depth = params$max_depth,
    learning_rate = params$eta,
    nrounds = params$nrounds,
    subsample = params$subsample,
    colsample_bytree = params$colsample_bytree,
    reg_lambda = params$lambda,        # ADD THIS
    reg_alpha = params$alpha,          # ADD THIS
    objective = "reg:squarederror",
  )  
  # Predictions on test sets
  pred1 <- predict(mynlm, newdata = as.matrix(X_matrix[sample_index1, , drop = FALSE]))
  pred2 <- predict(mynlm, newdata = as.matrix(X_matrix[sample_index2, , drop = FALSE]))
  
  y1 <- y[sample_index1]
  y2 <- y[sample_index2]
  
  # Calculate R² for both test sets
  R2orig1 <- calculate_r2(y1, pred1)
  R2orig2 <- calculate_r2(y2, pred2)
  R2_avg <- mean(c(R2orig1, R2orig2))
  
  # --- permutation-based drops ---
  Rnew1 <- sapply(seq_len(p), function(j)
    permR2(as.matrix(X_matrix[sample_index1, , drop = FALSE]), Y = y1, j = j, model = mynlm))
  Rnew2 <- sapply(seq_len(p), function(j)
    permR2(as.matrix(X_matrix[sample_index2, , drop = FALSE]), Y = y2, j = j, model = mynlm))
  
  beta1  <- R2orig1 - Rnew1
  beta2  <- R2orig2 - Rnew2
  mirror <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
  selected_index <- SelectFeatures(mirror, abs(mirror), q = 0.10)
  
  # Store results
  results <- rbind(results, data.frame(
    max_depth = params$max_depth,
    eta = params$eta,
    nrounds = params$nrounds,
    subsample = params$subsample,
    colsample_bytree = params$colsample_bytree,
    lambda = params$lambda,        # ADD THIS
    alpha = params$alpha,          # ADD THIS
    R2_test1 = R2orig1,
    R2_test2 = R2orig2,
    R2_avg = R2_avg,
    selected=length(selected_index)
  ))
  
  # Update best parameters
  if (selected > most_selected) {
    most_selected <- selected
    best_params <- params
    cat(sprintf("New best! max_depth=%d, eta=%.2f, nrounds=%d, subsample=%.1f, colsample_bytree=%.1f, R2_avg=%.4f\n",
                params$max_depth, params$eta, params$nrounds, params$subsample, params$colsample_bytree, R2_avg))
    print(selected)
  }
  
  # Progress update
  if (i %% 10 == 0) {
    cat(sprintf("Progress: %d/%d combinations tested\n", i, nrow(param_grid)))
  }
}

# Sort results by average R²
results <- results[order(-results$R2_avg), ]

# Display top 25 results
cat("\n=== Top 25 Parameter Combinations ===\n")
print(head(results, 25))

cat("\n=== Best Parameters ===\n")
cat(sprintf("max_depth: %d\n", best_params$max_depth))
cat(sprintf("eta: %.2f\n", best_params$eta))
cat(sprintf("nrounds: %d\n", best_params$nrounds))
cat(sprintf("subsample: %.1f\n", best_params$subsample))
cat(sprintf("colsample_bytree: %.1f\n", best_params$colsample_bytree))
cat(sprintf("Average R²: %.4f\n", best_r2))

# Train final model with best parameters
final_model <- xgboost(
  data = X_train,
  label = y_train,
  max_depth = best_params$max_depth,
  eta = best_params$eta,
  nrounds = best_params$nrounds,
  subsample = best_params$subsample,
  colsample_bytree = best_params$colsample_bytree,
  objective = "reg:squarederror",
  verbose = 0
)

# Final evaluation
pred1_final <- predict(final_model, newdata = as.matrix(X_matrix[sample_index1, , drop = FALSE]))
pred2_final <- predict(final_model, newdata = as.matrix(X_matrix[sample_index2, , drop = FALSE]))

R2_final1 <- calculate_r2(y[sample_index1], pred1_final)
R2_final2 <- calculate_r2(y[sample_index2], pred2_final)

cat("\n=== Final Model Performance (Best Parameters) ===\n")
cat(sprintf("R² on test set 1: %.4f\n", R2_final1))
cat(sprintf("R² on test set 2: %.4f\n", R2_final2))
cat(sprintf("Average R²: %.4f\n", mean(c(R2_final1, R2_final2))))

# Save results
write.csv(results, "hyperparameter_tuning_results.csv", row.names = FALSE)