library(randomForest)

# Function to calculate R²
calculate_r2 <- function(actual, predicted) {
  1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
}

# Define hyperparameter grid
# More comprehensive parameter grid
param_grid_extended <- expand.grid(
  degree = c(1, 2, 3, 4),
  nprune = seq(10, 80, by = 10),
  nk = c(50, 100, 200, 300),
  thresh = c(0.001, 0.01),       # Forward stepping threshold
  minspan = c(-1, 5, 10)         # Min span between knots (-1 = auto)
)

results_extended <- data.frame(
  degree = integer(),
  nprune = integer(),
  nk = integer(),
  thresh = numeric(),
  minspan = integer(),
  R2_test1 = numeric(),
  R2_test2 = numeric(),
  R2_avg = numeric(),
  n_terms = integer()
)

# prepare the data
data_full=mydata_full
X=mydata_full[,-c(1,2,3,4,8,10,11,12,13,18,19)]
#View(X)
y=mydata_full[,1]
# --- indices ---
set.seed(20260501)
amountTrain <- 0.5
amountTest  <- 1 - amountTrain
num_split   <- 1 #25#0#5010   # Number of data splits
q=0.1
n=nrow(X)
p=ncol(X)
train_index     <- sample.int(n, size = floor(0.5 * n), replace = FALSE)
remaining_index <- setdiff(seq_len(n), train_index)

# split the remaining part in two halves
size_half     <- floor((amountTest / 2) * n)
sample_index1 <- sample(remaining_index, size = size_half, replace = FALSE)
sample_index2 <- setdiff(remaining_index, sample_index1)
mydata_full=cbind(y,X)
dataTrain <- mydata_full[train_index, , drop = FALSE]
dataTrain
# Grid search
best_r2 <- -Inf
best_params <- NULL


for (i in 1:nrow(param_grid)) {
  params <- param_grid_extended[i, ]
  
  # Train model with current parameters
  mars_model <- earth(
    y ~ ., 
    data = dataTrain,
    degree = params$degree,
    nprune = params$nprune,
    nk = params$nk,
    thresh = params$thresh,
    minspan = params$minspan,
    pmethod = "backward"
  )
  
  pred1 <- predict(mars_model, newdata = as.data.frame(X[sample_index1, , drop = FALSE]))
  pred2 <- predict(mars_model, newdata = as.data.frame(X[sample_index2, , drop = FALSE]))
  
  
  y1 <- y[sample_index1]
  y2 <- y[sample_index2]
  
  # Calculate R² for both test sets
  R2orig1 <- calculate_r2(y1, pred1)
  R2orig2 <- calculate_r2(y2, pred2)
  R2_avg <- mean(c(R2orig1, R2orig2))
  if(R2_avg>0){
  
  # --- permutation-based drops ---
  Rnew1 <- sapply(seq_len(p), function(j)
    permR2(as.data.frame(X[sample_index1, , drop = FALSE]), Y = y1, j = j, model = mars_model))
  Rnew2 <- sapply(seq_len(p), function(j)
    permR2(as.data.frame(X[sample_index2, , drop = FALSE]), Y = y2, j = j, model = mars_model))
  
  beta1  <- R2orig1 - Rnew1
  beta2  <- R2orig2 - Rnew2
  mirror <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
  selected_index <- SelectFeatures(mirror, abs(mirror), q = 0.10)
  num_sel <- length(selected_index)
  }else{
    num_sel=-1
}
  # Store results
  results_extended <- rbind(results_extended, data.frame(
    degree = params$degree,
    nprune = params$nprune,
    nk = params$nk,
    thresh = params$thresh,
    minspan = params$minspan,
    R2_test1 = R2orig1,
    R2_test2 = R2orig2,
    R2_avg = R2_avg,
    n_terms =num_sel
  ))
  
  # Update best parameters
  if (R2_avg > best_r2) {
    best_r2 <- R2_avg
    best_params <- params
    cat(sprintf("New best! mtry=%d, ntree=%d, nodesize=%d, selected=%.4f\n",
                params$mtry, params$ntree, params$nodesize, num_sel))
  }
  
  # Progress update
  if (i %% 10 == 0) {
    cat(sprintf("Progress: %d/%d combinations tested\n", i, nrow(param_grid)))
  }
}

# Sort results by average R²
results_extended <- results_extended[order(-results_extended$selected), ]

# Display top 10 results
cat("\n=== Top 10 Parameter Combinations ===\n")
print(head(results, 25))

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