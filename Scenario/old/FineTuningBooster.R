library(xgboost)

# --- your fixed settings ------------------------
set.seed(456)
num_split <- 50
n <-500
p <- 100
p0 <- 10
q <- 0.1
signal_index <- sample(c(1:p), size = p0, replace = F)
delta=10
# simulate data
n1 <- floor(n/2); n2 <- n - n1
X1 <- matrix(rnorm(n1*p, mean= 1), n1, p)
X2 <- matrix(rnorm(n2*p, mean=-1), n2, p)
X  <- rbind(X1, X2)
beta_star <- numeric(p)
beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))*10
y <- (X^2 %*% beta_star + rnorm(n))

nrounds    <- 500
num_split  <- 1
amountTrain<- 0.333

# assume X (n×p) and y (length n) already exist
data <- data.frame(y = y, X)

# --- define a grid of candidate params -----------
param_grid <- expand.grid(
  eta               = c(0.01, 0.05, 0.1),
  max_depth         = c(3, 4, 6),
  subsample         = c(0.6, 0.8, 1.0),
  colsample_bytree  = c(0.6, 0.8, 1.0),
  lambda            = c(0,0.1, 0.5, 1, 5),
  alpha             = c(0,0.1, 0.5, 1),
  booster           = "gbtree",
  stringsAsFactors  = FALSE
)

# prepare storage
param_grid$mean_R2 <- NA_real_

# helper to compute R2
calc_r2 <- function(obs, pred) {
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
}

# --- grid search over params ----------------------
for(j in 1:nrow(param_grid)) {
  pars   <- as.list(param_grid[j, c("eta","max_depth","subsample",
                                    "colsample_bytree","lambda","alpha",
                                    "booster")])
  R2_vals <- numeric(num_split)
  
  for(iter in seq_len(num_split)) {
    # split train/test
    train_idx <- sample(seq_len(nrow(data)), size = amountTrain * nrow(data))
    test_idx  <- setdiff(seq_len(nrow(data)), train_idx)
    
    dtrain <- xgb.DMatrix(
      data  = as.matrix(data[train_idx, -1]),
      label = data$y[train_idx]
    )
    dtest  <- xgb.DMatrix(
      data  = as.matrix(data[test_idx,  -1]),
      label = data$y[test_idx]
    )
    
    # fit model
    bst <- xgb.train(
      params   = pars,
      data     = dtrain,
      nrounds  = nrounds,
      verbose  = 0
    )
    
    # predict & compute R2
    pred         <- predict(bst, dtest)
    R2_vals[iter] <- calc_r2(data$y[test_idx], pred)
  }
  
  # store average R2
  param_grid$mean_R2[j] <- mean(R2_vals)
  print(j)
}

# --- pick best -------------------------------
best_row   <- which.max(param_grid$mean_R2)
best_param <- param_grid[best_row, ]
print(best_param)

View(param_grid)
#eta max_depth subsample colsample_bytree lambda alpha booster   mean_R2
#38 417 0.1         3       0.8              0.6      0   0.1  gbtree 0.7306313