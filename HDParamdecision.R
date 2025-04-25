hyper_grid <- expand.grid(
  lambda = c(0, 1, 2, 5),
  alpha  = c(0, 0.1, 1, 2),
  eta    = c(0.01, 0.05, 0.1)
)

# Store CV results
search_results <- data.frame()

#-------------------------
# 2) Cross-validation to find the best combination
#-------------------------
for(i in seq_len(nrow(hyper_grid))) {
  
  params <- list(
    booster   = "gblinear",
    objective = "reg:squarederror",
    lambda    = hyper_grid$lambda[i],
    alpha     = hyper_grid$alpha[i],
    eta       = hyper_grid$eta[i]
  )
  
  # Perform 5-fold cross-validation
  cv_model <- xgb.cv(
    params               = params,
    data                 = data.matrix(Xtrain), 
    label                = y[train_index],
    nrounds              = 500,
    nfold                = 5,
    early_stopping_rounds = 20,  # helps avoid overfitting
    verbose              = FALSE
  )
  
  # Get the best iteration and RMSE from CV
  best_iter <- cv_model$best_iteration
  best_rmse <- cv_model$evaluation_log[best_iter, "test_rmse_mean"]
  
  # Track the results
  search_results <- rbind(
    search_results,
    data.frame(
      lambda       = params$lambda,
      alpha        = params$alpha,
      eta          = params$eta,
      best_iter    = best_iter,
      cv_test_rmse = best_rmse
    )
  )
}

#-------------------------
# 3) Identify the best hyperparameters
#-------------------------
best_params <- search_results[which.min(search_results$test_rmse_mean), ]
best_params
