# Inputs:
##   X: n x p numeric matrix/data.frame of predictors
##   y: length-n numeric response
## Packages needed: mgcv, e1071 (optional: energy for distance correlation)
library(e1071)
library(MASS )
library(corrplot )
set.seed(42)
n <- 40
p <- 50
p0 <- 25
q <- 0.1
delta=10
n1 <- floor(n/2); n2 <- n - n1
S <- toeplitz((p:1)/p)
X1 <- mvrnorm(n1*p,mu=rep(0,p),Sigma=S)
X2 <- mvrnorm(n2*p,mu=rep(0,p),Sigma=S)
X  <- rbind(X1, X2)
beta_star <- numeric(p)
signal_index <- sample(c(1:p), size = p0, replace = F)
beta_star[signal_index] <- 100#rnorm(p0, 0, delta*sqrt(log(p)/n))
y <- (X^2 %*% beta_star) # + rnorm(n))

cor(X)
# --- Split: 1/3 train, 1/3 val, 1/3 test --------------------------------
n <- nrow(X)
idx <- sample(seq_len(n))
n_train <- floor(n/3)
n_test  <- floor(n/3)
n_val   <- n - n_train - n_test

train_id <- idx[1:n_train]
val_id   <- idx[(n_train + 1):(n_train + n_val)]
test_id  <- idx[(n_train + n_val + 1):n]

# DMatrix
X_train <- as.matrix(X[train_id, , drop = FALSE])
X_val   <- as.matrix(X[val_id,   , drop = FALSE])
X_test  <- as.matrix(X[test_id,  , drop = FALSE])
y_train <- scale(y[train_id])
y_val <-  scale(y[val_id]); y_test <-  scale(y[test_id])



library(kernlab)
set.seed(123)
# Define a tuning grid
tune_grid <- expand.grid(
  degree=c(1,2,3),
  lambda = c(0.001, 0.01, 0.1, 1),
  scale  = c(0.5, 1, 2),
  offset = c(0, 1)
)

# Initialize results storage
results <- data.frame(degree = NA, lambda = NA, scale = NA, offset = NA, RMSE = NA)

best_rmse <- Inf
best_model <- NULL

# Grid search
for(i in 1:nrow(tune_grid)){
  params <- tune_grid[i, ]
  
  model <- ksvm(
    x = as.matrix(X_train),
    y = y_train,
    kernel = polydot(degree = params$degree, scale = params$scale, offset = params$offset),
    lambda = params$lambda
  )
  
  y_pred <- predict(model, as.matrix(X_val))  # validation set
  rmse <- sqrt(mean((y_val - y_pred)^2))
  
  results[i, ] <- c(params, rmse)
  
  if(rmse < best_rmse){
    best_rmse <- rmse
    best_model <- model
  }
}

# Best model
best_model
myparms=results[which.min(results$RMSE), ]
myparms
# Fit a polynomial kernel regression
krr_model <- ksvm(
  x = as.matrix(X_train),
  y = y_train,
  kernel = polydot(degree = myparms$degree, scale = myparms$scale, offset = myparms$offset), lambda = myparms$lambda
)

# Predictions
y_pred <- predict(krr_model, as.matrix(X_test))

# Evaluate performance
mse <- mean((y_test - y_pred)^2)
sst <- sum((y_test - mean(y_test))^2)
sse <- sum((y_test - y_pred)^2)
r2  <- 1 - sse/sst
cat("Test R^2:", r2, "\n")
