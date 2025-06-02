library(lightgbm)
library(MASS)
set.seed(1)
num_split <- 1
n <-1500
p <- 2000
p0 <- 200
q <- 0.1
delta <- 10

# simulate data
X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
beta_star <- numeric(p)
signal_index=1:p0
beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))*100
y <- scale(X^2 %*% beta_star + rnorm(n))

# --- your fixed settings ------------------------
nrounds     <- 500
num_split   <- 1
amountTrain <- 0.333

# assume X (n×p) and y (length n) already exist
data <- data.frame(y = y, X)

# --- define a grid of candidate params -----------
#param_grid <- expand.grid(
#  learning_rate     = c(0.01, 0.05, 0.1),
#  eta       = c(0.01,0.05,0.1,0.2),
  #max_depth         = c(3, 4, 6),
#  max_depth         = c(1,2,3),
#  subsample         = c(0.4,0.5,0.6, 0.8, 1.0),
#  feature_fraction  = c(0.6, 0.8, 1.0),
#  lambda_l2         = c(0, 1, 5),
#  lambda_l1         = c(0, 1),
#  stringsAsFactors  = FALSE
#)
#learning_rate  eta max_depth subsample feature_fraction lambda_l2 lambda_l1
#75           0.1 0.01         1       0.6              0.6         0         0
#mean_R2
#75 0.07516241

param_grid <- expand.grid(
  learning_rate     = c(0.1,0.2,0.25),
  eta       = c(0.001,0.01),
  #max_depth         = c(3, 4, 6),
  max_depth         = c(1,3,6),
  subsample         = c(0.6),
  feature_fraction  = c(0.1,0.25,0.6),
  lambda_l2         = c(0, 1, 5),
  lambda_l1         = c(0, 1),
  stringsAsFactors  = FALSE
)

# prepare storage
param_grid$mean_R2 <- NA_real_

# helper to compute R2
calc_r2 <- function(obs, pred) {
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
}

# --- grid search over params ----------------------
for(j in seq_len(nrow(param_grid))) {
  pars <- as.list(param_grid[j, ])
  pars$objective <- "regression"
  pars$metric    <- "l2"
  pars$verbose   <- -1
  
  R2_vals <- numeric(num_split)
  
  for(iter in seq_len(num_split)) {
    # split train/test
    train_idx <- sample(seq_len(nrow(data)), size = amountTrain * nrow(data))
    test_idx  <- setdiff(seq_len(nrow(data)), train_idx)
    
    dtrain <- lgb.Dataset(
      data  = as.matrix(data[train_idx, -1]),
      label = data$y[train_idx]
    )
    
    dtest <- as.matrix(data[test_idx, -1])
    
    # fit model
    model <- lgb.train(
      params   = pars,
      data     = dtrain,
      nrounds  = nrounds,
      verbose  = -1
    )
    
    # predict & compute R2
    pred           <- predict(model, dtest)
    R2_vals[iter]  <- calc_r2(data$y[test_idx], pred)
  }
  
  # store average R2
  param_grid$mean_R2[j] <- mean(R2_vals)
  print(c(j,mean(R2_vals)))
}

# --- pick best -------------------------------
best_row   <- which.max(param_grid$mean_R2)
best_param <- param_grid[best_row, ]
print(best_param)
learning_rate max_depth subsample feature_fraction lambda_l2 lambda_l1    mean_R2
349          0.01         6         1              0.6         1         1 -0.0218049
View(param_grid)
