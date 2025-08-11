# Inputs:
##   X: n x p numeric matrix/data.frame of predictors
##   y: length-n numeric response
## Packages needed: mgcv, e1071 (optional: energy for distance correlation)
library(BART)
set.seed(42)
n <-400
p <- 500
p0 <- 25
q <- 0.1
delta=10
n1 <- floor(n/2); n2 <- n - n1
X1 <- matrix(rnorm(n1*p, mean=0), n1, p)
X2 <- matrix(rnorm(n2*p, mean=0), n2, p)
X  <- rbind(X1, X2)
beta_star <- numeric(p)
signal_index <- sample(c(1:p), size = p0, replace = F)
beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))*100
y <- (X^2 %*% beta_star + rnorm(n))
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
y_test
# Standardize predictors using TRAIN stats (helps tree proposals with high-p)
scale_fit   <- function(A) list(mu = colMeans(A), sd = pmax(apply(A, 2, sd), 1e-8))
scale_apply <- function(A, s) sweep(sweep(A, 2, s$mu, "-"), 2, s$sd, "/")

sf <- scale_fit(X_train)
Xt <- scale_apply(X_train, sf)
Xv <- scale_apply(X_val,   sf)
Xs <- scale_apply(X_test,  sf)

# --- BART tuning ---------------------------------------------------------
suppressPackageStartupMessages(library(dbarts))
fit <- dbarts::bart(
  x.train = X_train, y.train = y_train,
  nskip   = 1000,    # burn-in
  ndpost  = 1000,    # posterior draws
  keeptrees = FALSE,
  verbose = FALSE
)
summary(y_train); sd(y_train)
anyNA(y_train); any(!is.finite(y_train))
colSums(!is.finite(X_train))[1:10]  # first few columns
sum(apply(X_train, 2, sd) == 0)     # zero-variance columns
wb <- wbart(x.train = Xt, y.train = y_train, x.test = Xv, ndpost=1000, nskip=1000, ntree=200)
pred_test <- predict(wb,newdata=X_test)
R2_test   <- R2(y_test, pred_test)

# Modest grids (keep runtime reasonable with n=133)
ntree_grid   <- c(50, 100, 200)       # number of trees
k_grid       <- c(1, 2, 3)            # prior shrinkage (larger = stronger shrinkage)
power_grid   <- c(2.0)                # tree depth prior beta (default ~2)
base_grid    <- c(0.95)               # tree depth prior alpha (default ~0.95)
sigdf_grid   <- c(3)                  # df for sigma prior
sigquant_grid<- c(0.90, 0.99)         # quantile for sigma prior (robust vs. diffuse)

best <- list(r2 = -Inf)

for (nt in ntree_grid) {
  for (k in k_grid) {
    for (pw in power_grid) {
      for (bs in base_grid) {
        for (sdof in sigdf_grid) {
          for (sq in sigquant_grid) {
            fit <- dbarts::bart(
              x.train = Xt, y.train = y_train,
              ntree   = nt,
              k       = k,
              power   = pw,
              base    = bs,
              sigdf   = sdof,
              sigquant= sq,
              nskip   = 1000,    # burn-in
              ndpost  = 1000,    # posterior draws
              keeptrees = FALSE,
              verbose = FALSE
            )
            pred_test <- predict(fit,newdata=X_test)
            r2_val   <- R2(y_test, pred_test)
            if (is.finite(r2_val) && r2_val > best$r2) {
              best <- list(r2 = r2_val, ntree = nt, k = k, power = pw, base = bs,
                           sigdf = sdof, sigquant = sq)
            }
          }
        }
      }
    }
  }
}

# --- Refit on Train+Val with best hyperparameters; evaluate on Test ------
Xtv  <- rbind(Xt, Xv)
y_tv <- c(y_train, y_val)

fit_final <- dbarts::bart(
  x.train = Xtv, y.train = y_tv,
  x.test  = Xs,
  ntree   = best$ntree,
  k       = best$k,
  power   = best$power,
  base    = best$base,
  sigdf   = best$sigdf,
  sigquant= best$sigquant,
  nskip   = 1500,
  ndpost  = 1500,
  keeptrees = FALSE,
  verbose = FALSE
)

pred_test <- as.vector(fit_final$yhat.test.mean)
R2_test   <- R2(y_test, pred_test)

cat(sprintf(
  "BART | ntree=%d k=%g power=%g base=%g sigdf=%d sigquant=%.2f | Val R2=%.3f | Test R2=%.3f\n",
  best$ntree, best$k, best$power, best$base, best$sigdf, best$sigquant, best$r2, R2_test
))