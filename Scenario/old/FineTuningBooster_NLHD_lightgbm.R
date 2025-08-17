# ---- packages ----
library(lightgbm)
library(doParallel)
library(foreach)

# --- your fixed settings ------------------------
set.seed(456)
n        <- 200
p        <- 240
p0       <- 25
delta    <- 10
n1       <- floor(n/2); n2 <- n - n1
X1       <- matrix(rnorm(n1*p, mean=-0.5), n1, p)
X2       <- matrix(rnorm(n2*p, mean=0.5), n2, p)
X        <- rbind(X1, X2)
signal_index <- sample(p, p0)
beta_star    <- numeric(p)
beta_star[signal_index] <- rnorm(p0, 0, delta * sqrt(log(p)/n)) * 10
y        <- scale(X^2 %*% beta_star + rnorm(n))
data     <- data.frame(y = y, X)

# CV settings
amountTrain     <- 0.333
num_split       <- 3

# --- LightGBM DART parameter grid (rough XGB -> LGB mapping) ---
# learning_rate ~ eta
# max_depth retained; num_leaves will be set to 2^max_depth inside the loop
# lambda_l1 ~ alpha; lambda_l2 ~ lambda
# feature_fraction ~ colsample_bytree; bagging_fraction ~ subsample; bagging_freq enables bagging
# DART-specific: boosting = "dart", drop_rate, skip_drop, xgboost_dart_mode
param_grid <- expand.grid(
  eta              = c(0.03,0.05,0.1,0.3),  # slow learning
  max_depth        = c(1,2, 3),               # shallow (prevents memorizing noise)
  min_child_weight = 10,          # larger leaves → more conservative
  colsample_bytree = 0.8,      # don’t go too low or you’ll miss key vars
  subsample        = 0.9,      # row subsampling for variance control
  gamma            = 3,            # split penalty (higher = fewer splits)
  lambda           = 1,           # L2: moderate
  alpha            = 25,      # L1: strong → favors sparsity
  tree_method      = "exact",               # tiny data
  objective        = "reg:squarederror",
  eval_metric      = "rmse",
  booster          = "gbtree",
  stringsAsFactors = FALSE
)
param_grid$mean_R2 <- NA_real_

# helper
calc_r2 <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)

# parallel backend
ncores <- 24
cl     <- makeCluster(ncores)
registerDoParallel(cl)

# Ensure LightGBM uses a sensible number of threads per worker j=1
# (optional; LightGBM is multi-threaded; you can also set "num_threads" in params)
lgb_threads <- max(1, floor(parallel::detectCores() / ncores))

# parallel loop
results <- foreach(j = seq_len(nrow(param_grid)), .combine = rbind,
                   .packages = "lightgbm") %dopar% {
                     
                     # Build params list for this config
                     pg   <- param_grid[j, ]
                     pars <- list(
                       eval_metric      = as.character(pg$eval_metric), # e.g., "rmse"
                       booster          = "gbtree",
                       eta              = pg$eta,                      # learning rate
                       max_depth        = pg$max_depth,
                       min_child_weight = pg$min_child_weight,
                       colsample_bytree = pg$colsample_bytree,
                       subsample        = pg$subsample,
                       gamma            = pg$gamma,
                       lambda           = pg$lambda,                   # L2 regularization
                       alpha            = pg$alpha,                    # L1 regularization
                       tree_method      = as.character(pg$tree_method)
               # threads per worker
                     )
                     R2_vals <- numeric(num_split)
                     
                     for (iter in seq_len(num_split)) {
                       # train/test split
                       idx       <- sample(nrow(data), size = amountTrain * nrow(data))
                       train_idx <- idx
                       test_idx  <- setdiff(seq_len(nrow(data)), idx)
                       
                       x_train <- data.matrix(data[train_idx, -1, drop = FALSE])
                       y_train <- data$y[train_idx]
                       x_test  <- data.matrix(data[test_idx,  -1, drop = FALSE])
                       y_test  <- data$y[test_idx]
                       
                       dtrain <- lgb.Dataset(data = x_train, label = y_train)
                       dvalid <- lgb.Dataset(data = x_test,  label = y_test)
                       
                       # train with early stopping (valid set)
                       bst <- lgb.train(
                         params              = pars,
                         data                = dtrain,
                          verbose             = -1
                       )
                       
                       pred <- predict(bst, x_test)
                       R2_vals[iter] <- calc_r2(y_test, pred)
                     }
                     
                     c(j = j, mean_R2 = mean(R2_vals))
                   }

# write results back to param_grid
param_grid$mean_R2[results[, "j"]] <- results[, "mean_R2"]

# clean up
stopCluster(cl)
registerDoSEQ()

# inspect top configs
head(param_grid[order(-param_grid$mean_R2), ], 5)
