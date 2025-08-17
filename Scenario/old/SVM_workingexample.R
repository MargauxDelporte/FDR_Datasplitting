# Inputs:
##   X: n x p numeric matrix/data.frame of predictors
##   y: length-n numeric response
## Packages needed: mgcv, e1071 (optional: energy for distance correlation)
library(e1071)
set.seed(42)
n <-2500
p <- 50
p0 <- 50
q <- 0.1
delta=10
n1 <- floor(n/2); n2 <- n - n1
X1 <- matrix(rnorm(n1*p, mean=0), n1, p)
X2 <- matrix(rnorm(n2*p, mean=0), n2, p)
X  <- rbind(X1, X2)
beta_star <- numeric(p)
signal_index <- sample(c(1:p), size = p0, replace = F)
beta_star[signal_index] <- 100#rnorm(p0, 0, delta*sqrt(log(p)/n))
y <- (X^2 %*% beta_star) # + rnorm(n))
# --- Split: 1/3 train, 1/3 val, 1/3 test -------------------------------

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
y_train <- y[train_id]; y_val <- y[val_id]; y_test <- y[test_id]


## --------------------------
## Standardize predictors
## --------------------------
scale_train <- scale(X_train)
center <- attr(scale_train, "scaled:center")
scale  <- attr(scale_train, "scaled:scale")

scale_fn <- function(X) {
  scale(X, center = center, scale = scale)
}

X_train <- scale_train
X_val   <- scale_fn(X_val)
X_test  <- scale_fn(X_test)

## --------------------------
## Hyperparameter tuning
## --------------------------
set.seed(42)
tuned <- tune.svm(
  x = X_train,
  y = y_train,
  kernel = "polynomial",
  cost   = 2^(-3:5),
  gamma  = 2^(-6:2),
  degree = 2,              # polynomial degree
  tunecontrol = tune.control(cross = 5)
)
cat("Best parameters:\n")
print(tuned$best.parameters)

## --------------------------
## Refit on train + val
## --------------------------
X_trainval <- rbind(X_train, X_val)
y_trainval <- c(y_train, y_val)

final_model <- svm(
  x = X_trainval,
  y = y_trainval,
  kernel = "polynomial",
  cost   = tuned$best.parameters$cost,
  gamma  = tuned$best.parameters$gamma,
  degree = 2
)
#fit <- svm(y ~ ., data = dat, scale = FALSE, kernel = "radial", cost = 5)

## --------------------------
## Evaluate on test set
## --------------------------
y_pred <- predict(final_model, X_test)

mse <- mean((y_test - y_pred)^2)
sst <- sum((y_test - mean(y_test))^2)
sse <- sum((y_test - y_pred)^2)
r2  <- 1 - sse/sst

cat("Test MSE:", mse, "\n")
cat("Test R^2:", r2, "\n")
