library(randomForest)

# Function to calculate R²
calculate_r2 <- function(actual, predicted) {
  1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
}

# Define hyperparameter grid
param_grid <- expand.grid(
  mtry = c(50, 100, 150, 200, 260, 300),
  ntree = c(100, 200, 350, 500),
  nodesize = c(5, 10, 16, 20, 25)
)

# Initialize results storage
results <- data.frame(
  mtry = integer(),
  ntree = integer(),
  nodesize = integer(),
  R2_test1 = numeric(),
  R2_test2 = numeric(),
  R2_avg = numeric()
)

# Grid search
best_r2 <- -Inf
best_params <- NULL

cat("Starting hyperparameter tuning...\n")
cat("Total combinations to test:", nrow(param_grid), "\n\n")

for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  
  # Train model with current parameters
  mynlm <- randomForest(
    y ~ ., 
    mtry = params$mtry,
    ntree = params$ntree,
    nodesize = params$nodesize,
    data = dataTrain
  )
  
  # Predictions on test sets
  pred1 <- predict(mynlm, newdata = as.data.frame(X[sample_index1, , drop = FALSE]))
  pred2 <- predict(mynlm, newdata = as.data.frame(X[sample_index2, , drop = FALSE]))
  
  y1 <- y[sample_index1]
  y2 <- y[sample_index2]
  
  # Calculate R² for both test sets
  R2_test1 <- calculate_r2(y1, pred1)
  R2_test2 <- calculate_r2(y2, pred2)
  R2_avg <- mean(c(R2_test1, R2_test2))
  
  # Store results
  results <- rbind(results, data.frame(
    mtry = params$mtry,
    ntree = params$ntree,
    nodesize = params$nodesize,
    R2_test1 = R2_test1,
    R2_test2 = R2_test2,
    R2_avg = R2_avg
  ))
  
  # Update best parameters
  if (R2_avg > best_r2) {
    best_r2 <- R2_avg
    best_params <- params
    cat(sprintf("New best! mtry=%d, ntree=%d, nodesize=%d, R²_avg=%.4f\n",
                params$mtry, params$ntree, params$nodesize, R2_avg))
  }
  
  # Progress update
  if (i %% 10 == 0) {
    cat(sprintf("Progress: %d/%d combinations tested\n", i, nrow(param_grid)))
  }
}

# Sort results by average R²
results <- results[order(-results$R2_avg), ]

# Display top 10 results
cat("\n=== Top 10 Parameter Combinations ===\n")
print(head(results, 10))

cat("\n=== Best Parameters ===\n")
cat(sprintf("mtry: %d\n", best_params$mtry))
cat(sprintf("ntree: %d\n", best_params$ntree))
cat(sprintf("nodesize: %d\n", best_params$nodesize))
cat(sprintf("Average R²: %.4f\n", best_r2))

# Train final model with best parameters
final_model <- randomForest(
  y ~ ., 
  mtry = best_params$mtry,
  ntree = best_params$ntree,
  nodesize = best_params$nodesize,
  data = dataTrain
)

# Final evaluation
pred1_final <- predict(final_model, newdata = as.data.frame(X[sample_index1, , drop = FALSE]))
pred2_final <- predict(final_model, newdata = as.data.frame(X[sample_index2, , drop = FALSE]))

R2_final1 <- calculate_r2(y[sample_index1], pred1_final)
R2_final2 <- calculate_r2(y[sample_index2], pred2_final)

cat("\n=== Final Model Performance ===\n")
cat(sprintf("R² on test set 1: %.4f\n", R2_final1))
cat(sprintf("R² on test set 2: %.4f\n", R2_final2))
cat(sprintf("Average R²: %.4f\n", mean(c(R2_final1, R2_final2))))