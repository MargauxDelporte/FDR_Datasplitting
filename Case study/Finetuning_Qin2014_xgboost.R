library(xgboost)

# Function to calculate R²
# Function to calculate R²
calculate_r2 <- function(actual, predicted) {
  1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
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
  train_rmse = numeric(),
  test_rmse = numeric(),
  test_r2 = numeric()
)
# --- check dataset ---
set.seed(20260501)
amountTrain <- 0.5
amountTest  <- 1 - amountTrain
num_split   <- 1 #25#0#5010   # Number of data splits
q=0.1
n=nrow(X)
train_index     <- sample.int(n, size = floor(0.5 * n), replace = FALSE)
remaining_index <- setdiff(seq_len(n), train_index)
dtrain <- xgb.DMatrix(data = as.matrix(X_matrix[train_index,]), label = y[train_index])
dtest <- xgb.DMatrix(data = as.matrix(X_matrix[remaining_index,]), label = y[remaining_index])
y_train=y[train_index]
y_test=y[remaining_index]
# Grid search
best_r2 <- -Inf
best_params <- NULL

cat("Starting XGBoost hyperparameter tuning...\n")
cat("Total combinations to test:", nrow(param_grid), "\n\n")

for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  
  # Set parameters
  param_list <- list(
    objective = "reg:squarederror",
    eta = params$eta,
    max_depth = params$max_depth,
    subsample = params$subsample,
    colsample_bytree = params$colsample_bytree,
    lambda = params$lambda,
    alpha = params$alpha
  )
  
  # Train model with early stopping
  set.seed(123)
  xgb_model <- xgb.train(
    params = param_list,
    data = dtrain,
    nrounds = params$nrounds,
    watchlist = list(train = dtrain, test = dtest),
    early_stopping_rounds = 20,
    verbose = 0
  )
  
  # Get predictions
  train_pred <- predict(xgb_model, dtrain)
  test_pred <- predict(xgb_model, dtest)
  
  # Calculate metrics
  train_rmse <- sqrt(mean((y_train - train_pred)^2))
  test_rmse <- sqrt(mean((y_test - test_pred)^2))
  test_r2 <- calculate_r2(y_test, test_pred)
  
  # Store results
  results <- rbind(results, data.frame(
    eta = params$eta,
    max_depth = params$max_depth,
    subsample = params$subsample,
    colsample_bytree = params$colsample_bytree,
    lambda = params$lambda,
    alpha = params$alpha,
    nrounds = params$nrounds,  # Use best iteration from early stopping
    train_rmse = train_rmse,
    test_rmse = test_rmse,
    test_r2 = test_r2
  ))
  
  # Update best parameters
  if (test_r2 > best_r2) {
    best_r2 <- test_r2
    best_params <- params
    best_params$nrounds <- xgb_model$best_iteration
    cat(sprintf("New best! eta=%.3f, max_depth=%d, R²=%.4f, nrounds=%d\n",
                params$eta, params$max_depth, test_r2, xgb_model$best_iteration))
  }
  
  # Progress update
  if (i %% 10 == 0) {
    cat(sprintf("Progress: %d/%d combinations tested\n", i, nrow(param_grid)))
  }
}

# Sort results by test R²
results <- results[order(-results$test_r2), ]

# Display top 10 results
cat("\n=== Top 10 Parameter Combinations ===\n")
print(head(results, 10))

cat("\n=== Best Parameters ===\n")
print(best_params)
cat(sprintf("Best Test R²: %.4f\n", best_r2))
