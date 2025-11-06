# --- Strongly nonlinear scenario: XGBoost >> LASSO ---
# install.packages(c("xgboost","glmnet"))  # if needed
### High dimension linear model
rm(list = ls())

mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
#mywd='C:/Users/mde4023/Documents/GitHub/FDR_Datasplitting'

setwd(mywd)

source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Functions/TriangleBoosterTrainMS3.R'))
source(paste0(mywd,'/Functions/ApplyGBMKnockoff.R'))
#source(paste0(mywd,'/Functions/TriangleGBMTrainMS.R'))

source(paste0(mywd,'/Functions Dai/knockoff.R'))
source(paste0(mywd,'/Functions Dai/analysis.R'))
source(paste0(mywd,'/Functions Dai/MBHq.R'))
source(paste0(mywd,'/Functions Dai/DS.R'))
source(paste0(mywd,'/Functions Dai/fdp_power.R'))
library(xgboost)
library(glmnet)

set.seed(2025)

n <- 800
p <- 1000
p_true <- 20

# Simulate features
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("x", 1:p)

relu <- function(z) pmax(0, z)

term_step        <- function(x) as.numeric(x >  0.30)      # threshold step
term_relu_right  <- function(x) relu(x - 0.50)             # right ReLU
term_relu_left   <- function(x) relu(-x - 0.40)            # left ReLU
term_quad_ind    <- function(x) as.numeric((x^2) > 1.00)   # quadratic indicator
term_sin         <- function(x) sin(3.0 * x)               # smooth sine

# Build a signal from the first 20 covariates using ONE chosen transform
build_signal <- function(X, f, cols = 1:20, standardize = TRUE) {
  Z_list <- vector("list", length(cols))
  for (k in seq_along(cols)) {
    z <- f(X[, cols[k]])
    if (standardize) z <- scale(z)[, 1]   # unit variance per term
    Z_list[[k]] <- z
  }
  Reduce(`+`, Z_list)
}

# Common effect size and shared noise (so scenarios are directly comparable)
beta  <- 2.0
n     <- nrow(X)
noise <- rnorm(n, sd = 0.5)

# Signals (each uses ALL first 20 covariates with the SAME effect size)
signal_step       <- build_signal(X, term_step)
signal_relu_right <- build_signal(X, term_relu_right)
signal_relu_left  <- build_signal(X, term_relu_left)
signal_quad_ind   <- build_signal(X, term_quad_ind)
signal_sin        <- build_signal(X, term_sin)

# Responses (reuse the same noise for fair comparison)
y_step       <- beta * signal_step       + noise
y_relu_right <- beta * signal_relu_right + noise
y_relu_left  <- beta * signal_relu_left  + noise
y_quad_ind   <- beta * signal_quad_ind   + noise
y_sin        <- beta * signal_sin        + noise

# (Rest of your pipeline: split data, fit XGBoost vs LASSO, etc.)
signal_index=1:20
params <- list(
  objective = "reg:squarederror",
  max_depth = 4,          # allow higher-order interactions
  min_child_weight = 1,
  eta = 0.05,             # smaller step size
  subsample = 0.8,
  colsample_bytree = 0.8,
  reg_lambda = 1.0,
  reg_alpha = 0.0,
  nthread = 0             # use all CPU threads
)

g1 <- ApplyTriangleBoostTrain2( X = X, y = y_step, q = 0.1, num_split = 3,
                                signal_index = signal_index, myseed = 1)
g2 <- ApplyTriangleBoostTrain2( X = X, y = y_relu_right, q = 0.1, num_split = 3,
                                signal_index = signal_index, myseed = 1)
g3 <- ApplyTriangleBoostTrain2( X = X, y = y_relu_left, q = 0.1, num_split = 3,
                                signal_index = signal_index, myseed = 1)
g4 <- ApplyTriangleBoostTrain2( X = X, y = y_quad_ind, q = 0.1, num_split = 3,
                                signal_index = signal_index, myseed = 1)
g5 <- ApplyTriangleBoostTrain2( X = X, y = y_sin, q = 0.1, num_split = 3,
                                signal_index = signal_index, myseed = 1)

# (Rest of your pipeline: split data, fit XGBoost vs LASSO, etc.)

# Train/Valid/Test split: 60/20/20
n_train <- floor(0.6 * n)
n_valid <- floor(0.2 * n)
idx <- sample(1:n, n)
id_tr <- idx[1:n_train]
id_va <- idx[(n_train + 1):(n_train + n_valid)]
id_te <- idx[(n_train + n_valid + 1):n]

Xtr <- X[id_tr, , drop = FALSE]; ytr <- y[id_tr]
Xva <- X[id_va, , drop = FALSE]; yva <- y[id_va]
Xte <- X[id_te, , drop = FALSE]; yte <- y[id_te]

# ------------------- XGBoost -------------------
dtr <- xgb.DMatrix(Xtr, label = ytr)
dva <- xgb.DMatrix(Xva, label = yva)
dte <- xgb.DMatrix(Xte, label = yte)

params <- list(
  objective = "reg:squarederror",
  max_depth = 4,          # allow higher-order interactions
  min_child_weight = 1,
  eta = 0.05,             # smaller step size
  subsample = 0.8,
  colsample_bytree = 0.8,
  reg_lambda = 1.0,
  reg_alpha = 0.0,
  nthread = 0             # use all CPU threads
)

watchlist <- list(train = dtr, valid = dva)
xgb_model <- xgb.train(
  params = params,
  data = dtr,
  nrounds = 2000,                 # generous cap
  watchlist = watchlist,
  early_stopping_rounds = 50,     # stop when valid loss plateaus
  verbose = 0
)

xgb_pred <- predict(xgb_model, dte)
xgb_rmse <- sqrt(mean((xgb_pred - yte)^2))

# ------------------- LASSO -------------------
lasso_cv <- cv.glmnet(Xtr, ytr, alpha = 1, nfolds = 5, standardize = TRUE)
lasso_pred <- as.numeric(predict(lasso_cv, newx = Xte, s = "lambda.min"))
lasso_rmse <- sqrt(mean((lasso_pred - yte)^2))

cat("Strongly nonlinear scenario RMSEs:\n")
cat(sprintf("  XGBoost: %.3f\n", xgb_rmse))
cat(sprintf("  LASSO  : %.3f\n", lasso_rmse))

# Optional: show top features XGBoost used
#imp <- xgb.importance(model = xgb_model, feature_names = colnames(X))
#print(head(imp, 12))