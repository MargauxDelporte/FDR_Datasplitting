library(xgboost)

# --- your fixed settings ------------------------
set.seed(456)
num_split <- 5
n <-400
p <- 500
p0 <- 25
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
#param_grid <- expand.grid(
#  eta               = c(0.01, 0.05, 0.1,0.3),
#  max_depth         = c(4,5,6,8),
#  subsample         = c(0.6, 0.8, 1.0),
#  colsample_bytree  = c(0.6, 0.8, 1.0),
# lambda            = c(0,0.1, 0.5, 1, 5),
#  alpha             = c(0,0.1, 0.5, 1,1.5),
##  booster           = "gbtree",
#  stringsAsFactors  = FALSE
#)
param_grid <- expand.grid(
  eta               = c(0.005, 0.01, 0.05, 0.1),
  max_depth         = c(6, 8, 10, 12, 15),
  min_child_weight  = c(1, 5, 10),
  gamma             = c(0, 0.1, 1, 5),
  subsample         = c(0.4, 0.6, 0.8, 1),
  colsample_bytree  = c(0.4, 0.6, 0.8, 1),
  lambda            = c(0, 0.1, 1),
  alpha             = c(0, 0.1, 1),
  stringsAsFactors  = FALSE
)
# prepare storage
param_grid$mean_R2 <- NA_real_

# helper to compute R2
calc_r2 <- function(obs, pred) {
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
}

# --- grid search over params ---------------------- j=1
for(j in 1:nrow(param_grid)) {
  pars   <- as.list(param_grid[j, c("eta","max_depth",'min_child_weight','gamma',
                                    "subsample",
                                    "colsample_bytree","lambda","alpha")])
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
      params              = pars,
      data                = dtrain,
      nrounds             = 2000,                # high upper bound
      early_stopping_rounds = 50,                # stops if no improvement
      watchlist           = list(train = dtrain, eval = dtest),
      metric              = "rmse",              
      verbose             = 0
    )
    
    # predict & compute R2
    pred         <- predict(bst, dtest)
    R2_vals[iter] <- calc_r2(data$y[test_idx], pred)
  }
  
  # store average R2
  param_grid$mean_R2[j] <- mean(R2_vals)
  print(j)
}

cl <- makeCluster(20)
results_list <- foreach(
  k = seq_len(nrow(param_grid)),
  .packages = pkgs,
  .combine  = rbind
) %dopar% {
  s_val <- param_grid$s[k]

  # compute chunk of results
  chunk <- param_grid(i = i_val, s = s_val)
  
  # write out this chunk immediately
  fname <- sprintf("Results_s%02d_i%02d.csv", s_val, i_val)
  write.csv(chunk, file = paste0(mywd,"/Temp/",fname), row.names = FALSE)
  
  # return for final binding
  chunk
}
# --- pick best -------------------------------
best_row   <- which.max(param_grid$mean_R2)
best_param <- param_grid[best_row, ]
print(best_param)

View(param_grid)
#eta max_depth subsample colsample_bytree lambda alpha booster   mean_R2
#668 0.05         3         1              0.6      1   0.1  gbtree 0.4305495 if n=100 p=150
#20 0.05         3         1              0.6      0     0  gbtree 0.3945439 if n=200 and p=250
# 111 0.1         4       0.8              0.8    0.1     0  gbtree 0.424740 if n=100 p=150
#1182 0.05         4       0.6              0.8    0.1     1  gbtree 0.3047055 if n=400 and p=500