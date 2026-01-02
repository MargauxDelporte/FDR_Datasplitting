library(earth)
library(caret)

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

best_r2 <- -Inf
best_params <- NULL

cat("Starting extended MARS tuning...\n")

for (i in 1:min(100, nrow(param_grid_extended))) {  # Limit to 100 iterations
  params <- param_grid_extended[i, ]
  
  tryCatch({
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
    
    R2_test1 <- calculate_r2(y[sample_index1], pred1)
    R2_test2 <- calculate_r2(y[sample_index2], pred2)
    R2_avg <- mean(c(R2_test1, R2_test2))
    
    results_extended <- rbind(results_extended, data.frame(
      degree = params$degree,
      nprune = params$nprune,
      nk = params$nk,
      thresh = params$thresh,
      minspan = params$minspan,
      R2_test1 = R2_test1,
      R2_test2 = R2_test2,
      R2_avg = R2_avg,
      n_terms = length(mars_model$selected.terms)
    ))
    
    if (R2_avg > best_r2) {
      best_r2 <- R2_avg
      best_params <- params
      cat(sprintf("Iteration %d: New best! R²_avg=%.4f (terms=%d)\n", 
                  i, R2_avg, length(mars_model$selected.terms)))
    }
    
  }, error = function(e) {
    # Skip failed combinations
  })
}

results_extended <- results_extended[order(-results_extended$R2_avg), ]
print(head(results_extended, 15))