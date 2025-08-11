library(xgboost)
library(doParallel)
library(foreach)

# --- your fixed settings ------------------------
set.seed(456)
n        <- 500
p        <- 400
p0       <- 25
delta    <- 10
n1       <- floor(n/2); n2 <- n - n1
X1       <- matrix(rnorm(n1*p, mean= 0.5), n1, p)
X2       <- matrix(rnorm(n2*p, mean=-0.5), n2, p)
X        <- rbind(X1, X2)
signal_index <- sample(p, p0)
beta_star    <- numeric(p)
beta_star[signal_index] <- rnorm(p0, 0, delta * sqrt(log(p)/n)) * 10
y        <- (X^2 %*% beta_star + rnorm(n))
data     <- data.frame(y = y, X)

# CV settings
nrounds      <- 2000
early_stop   <- 50
amountTrain  <- 0.333
num_split    <- 3

# param grid
param_grid <- expand.grid(
  eta               = c(0.03, 0.10),   # learning rate
  max_depth         = c(4, 8),
  min_child_weight  = c(1, 5),
  subsample         = 0.9,              # fix to keep grid small
  colsample_bytree  = c(0.8, 1.0),
  alpha             = c(0, 1),         # L1
  lambda            = 1,                # L2
  rate_drop         = c(0.05, 0.15),   # DART dropout rate
  skip_drop         = 0.5,              # prob to skip dropout entirely
  sample_type       = c("uniform","weighted"),
  normalize_type    = "tree",           # "forest" also possible
  one_drop          = 0,                # try 1 in a follow-up sweep if desired
  tree_method       = "hist",           # "exact" if dataset is small
  objective         = "reg:squarederror",
  eval_metric       = "rmse",
  booster           = "dart",
  stringsAsFactors  = FALSE
)
param_grid$mean_R2 <- NA_real_

# helper
calc_r2 <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)

# set up parallel backend
ncores <- 24
cl     <- makeCluster(ncores)
registerDoParallel(cl)

# parallel loop
results <- foreach(j = seq_len(nrow(param_grid)), .combine = rbind, 
                   .packages = "xgboost") %dopar% {
                     
                     pars   <- as.list(param_grid[j, c("eta","max_depth",
                                                       #"min_child_weight",
                                                       #"gamma","subsample","colsample_bytree",
                                                       "lambda","alpha")])
                     R2_vals <- numeric(num_split)
                     
                     for (iter in seq_len(num_split)) {
                       # train/test split
                       idx       <- sample(nrow(data), size = amountTrain * nrow(data))
                       train_idx <- idx
                       test_idx  <- setdiff(seq_len(nrow(data)), idx)
                       
                       dtrain <- xgb.DMatrix(data = as.matrix(data[train_idx, -1]),
                                             label = data$y[train_idx])
                       dtest  <- xgb.DMatrix(data = as.matrix(data[test_idx,  -1]),
                                             label = data$y[test_idx])
                       
                       # fit with early stopping
                       bst <- xgb.train(
                         params                = pars,
                         data                  = dtrain,
                         nrounds               = nrounds,
                         early_stopping_rounds = early_stop,
                         watchlist             = list(eval = dtest),
                         metric                = "rmse",
                         verbose               = 0
                       )
                       
                       pred       <- predict(bst, dtest)
                       R2_vals[iter] <- calc_r2(data$y[test_idx], pred)
                     }
                     
                     # return j and mean R2
                     c(j = j, mean_R2 = mean(R2_vals))
                   }

# write results back to param_grid
param_grid$mean_R2[results[,"j"]] <- results[,"mean_R2"]

# clean up
stopCluster(cl)
registerDoSEQ()

# inspect top configs
head(param_grid[order(-param_grid$mean_R2), ], 5)





