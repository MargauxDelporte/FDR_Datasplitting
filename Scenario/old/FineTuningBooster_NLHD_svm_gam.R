# Inputs:
##   X: n x p numeric matrix/data.frame of predictors
##   y: length-n numeric response
## Packages needed: mgcv, e1071 (optional: energy for distance correlation)

set.seed(42)
n <-800
p <- 1000
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
n <- nrow(X)
idx <- sample(seq_len(n))
n_train <- floor(n/3)
n_test  <- floor(n/3)
n_val   <- n - n_train - n_test

train_id <- idx[1:n_train]
val_id   <- idx[(n_train + 1):(n_train + n_val)]
test_id  <- idx[(n_train + n_val + 1):n]

X_train <- as.data.frame(X[train_id, , drop = FALSE])
X_val   <- as.data.frame(X[val_id,   , drop = FALSE])
X_test  <- as.data.frame(X[test_id,  , drop = FALSE])
y_train <- y[train_id]; y_val <- y[val_id]; y_test <- y[test_id]

# Consistent column names
names(X_train) <- names(X_val) <- names(X_test) <- paste0("V", seq_len(ncol(X_train)))

# Helper: R^2
R2 <- function(y_true, y_pred) 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2)

# ======================================================================
# 1) GAM (mgcv) with shrinkage smooths; tuning k_basis and #features
# ======================================================================
library(mgcv)

# Nonlinear screening on TRAIN ONLY
dcor_screen <- function(Xdf, y) {
  if (requireNamespace("energy", quietly = TRUE)) {
    sapply(Xdf, function(z) energy::dcor(z, y))
  } else {
    # Fallback: correlation with squared predictor (captures quadratic signal)
    sapply(Xdf, function(z) abs(suppressWarnings(cor(z^2, y))))
  }
}
scores  <- dcor_screen(X_train, y_train)
ranked  <- order(scores, decreasing = TRUE)

k_basis_grid <- c(3, 4)
keep_grid    <- c(10, 20, 30)  # how many predictors to keep
best_gam <- list(score = -Inf)

for (k_basis in k_basis_grid) {
  for (k_keep in keep_grid) {
    keep_idx <- ranked[seq_len(min(k_keep, ncol(X_train)))]
    Xtr <- X_train[, keep_idx, drop = FALSE]
    Xv  <- X_val[,   keep_idx, drop = FALSE]
    
    dat_tr <- data.frame(y = y_train, Xtr)
    terms  <- paste(sprintf('s(%s, k=%d, bs="ts")', names(Xtr), k_basis), collapse = " + ")
    form   <- as.formula(paste("y ~", terms))
    
    fit <- gam(form, data = dat_tr, method = "REML", select = TRUE)
    pred_val <- predict(fit, newdata = Xv, type = "response")
    r2_val <- R2(y_val, pred_val)
    
    if (is.finite(r2_val) && r2_val > best_gam$score) {
      best_gam <- list(score = r2_val, k_basis = k_basis, k_keep = k_keep,
                       keep_idx = keep_idx, fit = fit)
    }
  }
}

# Refit GAM on train+val with best settings; evaluate on test
keep_idx <- best_gam$keep_idx
Xtrv     <- rbind(X_train[, keep_idx, drop = FALSE], X_val[, keep_idx, drop = FALSE])
y_trv    <- c(y_train, y_val)
dat_trv  <- data.frame(y = y_trv, Xtrv)
terms    <- paste(sprintf('s(%s, k=%d, bs="ts")', names(Xtrv), best_gam$k_basis), collapse = " + ")
form     <- as.formula(paste("y ~", terms))

fit_gam_final <- gam(form, data = dat_trv, method = "REML", select = TRUE)
pred_test_gam <- predict(fit_gam_final, newdata = X_test[, keep_idx, drop = FALSE])
R2_test_gam   <- R2(y_test, pred_test_gam)

# ======================================================================
# 2) SVR with degree-2 polynomial kernel; tuning hyperparameters
# ======================================================================
library(e1071)

gamma_grid <- c(0.5 / ncol(X_train), 1 / ncol(X_train), 2 / ncol(X_train))
coef0_grid <- c(0, 1)
cost_grid  <- c(0.1, 1, 10, 100)
eps_grid   <- c(0.01, 0.1, 0.3)

best_svr <- list(score = -Inf)

for (g in gamma_grid) {
  for (c0 in coef0_grid) {
    for (C in cost_grid) {
      for (eps in eps_grid) {
        fit <- svm(x = X_train, y = y_train,
                   type = "eps-regression",
                   kernel = "polynomial",
                   degree = 2,
                   gamma = g, coef0 = c0,
                   cost = C, epsilon = eps,
                   scale = TRUE)
        pred_val <- predict(fit, X_val)
        r2_val <- R2(y_val, pred_val)
        
        if (is.finite(r2_val) && r2_val > best_svr$score) {
          best_svr <- list(score = r2_val, gamma = g, coef0 = c0,
                           cost = C, epsilon = eps, fit = fit)
        }
      }
    }
  }
}

# Refit SVR on train+val with best hyperparams; evaluate on test
svr_final <- svm(x = rbind(X_train, X_val),
                 y = c(y_train, y_val),
                 type   = "eps-regression",
                 kernel = "polynomial",
                 degree = 2,
                 gamma  = best_svr$gamma,
                 coef0  = best_svr$coef0,
                 cost   = best_svr$cost,
                 epsilon= best_svr$epsilon,
                 scale  = TRUE)

pred_test_svr <- predict(svr_final, X_test)
R2_test_svr   <- R2(y_test, pred_test_svr)

# ---------------------- Results ---------------------------------------
cat(sprintf("GAM  | best: k_basis=%d, kept_vars=%d | Val R2=%.3f | Test R2=%.3f\n",
            best_gam$k_basis, best_gam$k_keep, best_gam$score, R2_test_gam))
cat(sprintf("SVR2 | best: gamma=%.4g, coef0=%g, cost=%g, epsilon=%g | Val R2=%.3f | Test R2=%.3f\n",
            best_svr$gamma, best_svr$coef0, best_svr$cost, best_svr$epsilon, best_svr$score, R2_test_svr))
