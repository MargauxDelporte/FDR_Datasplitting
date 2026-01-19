library(xgboost)


X_matrix=read.csv('X_matrix_Qin2014.csv')
y=read.csv('y.csv')

# Function to calculate R²
calculate_r2 <- function(actual, predicted) {
  1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
}

# Define hyperparameter grid for XGBoost
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
  train_rmse = numeric(),
  test_rmse = numeric(),
  test_r2 = numeric(),
  selected = numeric()
)

names(mydata)
# prepare the data
nrow(mydata_full) #237

# --- indices ---
set.seed(20260501)
amountTrain <- 0.5
amountTest  <- 1 - amountTrain
num_split   <- 1 #25#0#5010   # Number of data splits
q <- 0.1
n <- nrow(X)
train_index     <- sample.int(n, size = floor(0.5 * n), replace = FALSE)
remaining_index <- setdiff(seq_len(n), train_index)

# split the remaining part in two halves
size_half     <- floor((amountTest / 2) * n)
sample_index1 <- sample(remaining_index, size = size_half, replace = FALSE)
sample_index2 <- setdiff(remaining_index, sample_index1)

# Prepare data for XGBoost (requires matrix format)
X_train <- X_matrix[train_index,]
y_train <- y[train_index]

# Grid search
best_r2 <- -Inf
best_params <- NULL

cat("Starting hyperparameter tuning...\n")
cat("Total combinations to test:", nrow(param_grid), "\n\n")

for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  
  # Train model with current parameters
  mynlm <- xgboost(
    data = X_train,
    label = y_train,
    max_depth = params$max_depth,
    eta = params$eta,
    nrounds = params$nrounds,
    subsample = params$subsample,
    colsample_bytree = params$colsample_bytree,
    objective = "reg:squarederror",
    verbose = 0
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
    permR2(as.matrix(X[sample_index1, , drop = FALSE]), Y = y1, j = j, model = mynlm))
  Rnew2 <- sapply(seq_len(p), function(j)
    permR2(as.matrix(X[sample_index2, , drop = FALSE]), Y = y2, j = j, model = mynlm))
  
  beta1  <- R2orig1 - Rnew1
  beta2  <- R2orig2 - Rnew2
  mirror <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
  selected_index <- SelectFeatures(mirror, abs(mirror), q = 0.10)
  num_sel <- length(selected_index)
  
  # Store results
  results <- rbind(results, data.frame(
    max_depth = params$max_depth,
    eta = params$eta,
    nrounds = params$nrounds,
    subsample = params$subsample,
    colsample_bytree = params$colsample_bytree,
    R2_test1 = R2orig1,
    R2_test2 = R2orig2,
    R2_avg = R2_avg,
    selected = num_sel
  ))
  
  # Update best parameters
  if (R2_avg > best_r2) {
    best_r2 <- R2_avg
    best_params <- params
    cat(sprintf("New best! max_depth=%d, eta=%.2f, nrounds=%d, subsample=%.1f, colsample_bytree=%.1f, selected=%d\n",
                params$max_depth, params$eta, params$nrounds, params$subsample, params$colsample_bytree, num_sel))
  }
  
  # Progress update
  if (i %% 10 == 0) {
    cat(sprintf("Progress: %d/%d combinations tested\n", i, nrow(param_grid)))
  }
}

# Sort results by average R²
results <- results[order(-results$selected), ]
results <- results[order(-results$R2_avg), ]

# Display top 10 results
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
pred1_final <- predict(final_model, newdata = as.matrix(X[sample_index1, , drop = FALSE]))
pred2_final <- predict(final_model, newdata = as.matrix(X[sample_index2, , drop = FALSE]))

R2_final1 <- calculate_r2(y[sample_index1], pred1_final)
R2_final2 <- calculate_r2(y[sample_index2], pred2_final)

cat("\n=== Final Model Performance (Best Parameters) ===\n")
cat(sprintf("R² on test set 1: %.4f\n", R2_final1))
cat(sprintf("R² on test set 2: %.4f\n", R2_final2))
cat(sprintf("Average R²: %.4f\n", mean(c(R2_final1, R2_final2))))

# Train final model with manually specified parameters (if needed)
final_model_manual <- xgboost(
  data = X_train,
  label = y_train,
  max_depth = 5,
  eta = 0.1,
  nrounds = 500,
  subsample = 0.8,
  colsample_bytree = 0.8,
  objective = "reg:squarederror",
  verbose = 0
)

# Final evaluation with manual parameters
pred1_manual <- predict(final_model_manual, newdata = as.matrix(X[sample_index1, , drop = FALSE]))
pred2_manual <- predict(final_model_manual, newdata = as.matrix(X[sample_index2, , drop = FALSE]))

R2_manual1 <- calculate_r2(y[sample_index1], pred1_manual)
R2_manual2 <- calculate_r2(y[sample_index2], pred2_manual)

cat("\n=== Final Model Performance (Manual Parameters) ===\n")
cat(sprintf("R² on test set 1: %.4f\n", R2_manual1))
cat(sprintf("R² on test set 2: %.4f\n", R2_manual2))
cat(sprintf("Average R²: %.4f\n", mean(c(R2_manual1, R2_manual2))))



