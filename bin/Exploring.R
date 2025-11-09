# ================================================
# run_xgb_fixed_params.R
# - Strongly nonlinear scenario (step signal)
# - Uses fixed XGBoost params picked from prior grid search
# - Prints test RMSEs and top-20 feature importances
# ================================================

suppressPackageStartupMessages({
  library(xgboost)
  library(glmnet)
})

set.seed(2025)

# -------------------- Data --------------------
n <- 800
p <- 1000
p_true <- 20
beta <- 1
signal_index=1:p_true
# Simulate features
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("x", 1:p)

# Step signal on first 20 covariates i=0
result=rep(0,5)
for(i in seq(from=0,to=2.9,by=0.1)){
S <- scale((X[, 1:p_true] > i) * 1)    # 0/1 -> numeric -> standardize per column
signal_step <- rowSums(S)
y <- beta *signal_step + rnorm(n)

g1 <- ApplyTriangleBoostTrain2( X = X, y = y, q = q, num_split = 2,signal_index = signal_index, myseed = 1)
DS_result      <- DS(          X = X, y = y, q = q, num_split = 2)

result=rbind(result,c(i,g1$DS_power,g1$DS_fdp,DS_result$DS_power,DS_result$DS_fdp))
print(i)
}
result=as.data.frame(result)
names(result)=c('seed','mypower','myfdp','dspower','dsfdp')
# -------------------- Split --------------------
# Train/Valid/Test = 60/20/20
n_train <- floor(0.6 * n)
n_valid <- floor(0.2 * n)

idx <- sample.int(n, n)
id_tr <- idx[1:n_train]
id_va <- idx[(n_train + 1):(n_train + n_valid)]
id_te <- idx[(n_train + n_valid + 1):n]

Xtr <- X[id_tr, , drop = FALSE]; ytr <- y[id_tr]
Xva <- X[id_va, , drop = FALSE]; yva <- y[id_va]
Xte <- X[id_te, , drop = FALSE]; yte <- y[id_te]

# -------------------- XGBoost (fixed best params) --------------------
# From your grid result:
# max_depth=2, min_child_weight=10, eta=0.2, subsample=0.6,
# colsample_bytree=0.6, reg_lambda=0, reg_alpha=1, best_nrounds=88
best_nrounds <- 88L

dtrva <- xgb.DMatrix(rbind(Xtr, Xva), label = c(ytr, yva))
dte   <- xgb.DMatrix(Xte, label = yte)

params_best <- list(
  objective        = "reg:squarederror",
  eval_metric      = "rmse",
  max_depth        = 2L,
  min_child_weight = 10,
  eta              = 0.2,
  subsample        = 0.6,
  colsample_bytree = 0.6,
  reg_lambda       = 0,
  reg_alpha        = 1,
  # max_depth            = 2,
  # eta                  = 0.1,
  # gamma                = 0,
  # colsample_bytree     = 0.6,
  # subsample            = 0.8,
  # min_child_weight     = 5,
  # reg_lambda           = 10,
  # reg_alpha            = 0,
  nthread          = 0  # all cores
)

fit_best <- xgb.train(
  params  = params_best,
  data    = dtrva,
  nrounds = best_nrounds,
  verbose = 0
)

pred_xgb <- predict(fit_best, dte)
rmse_xgb <- sqrt(mean((pred_xgb - yte)^2))

# -------------------- LASSO baseline --------------------
lasso_cv <- cv.glmnet(Xtr, ytr, alpha = 1, nfolds = 5, standardize = TRUE)
pred_lasso <- as.numeric(predict(lasso_cv, newx = Xte, s = "lambda.min"))
rmse_lasso <- sqrt(mean((pred_lasso - yte)^2))

cat("======== Test RMSE ========\n")
cat(sprintf("XGBoost (fixed params): %.4f\n", rmse_xgb))
cat(sprintf("LASSO               : %.4f\n\n", rmse_lasso))

# -------------------- Top 20 features --------------------
imp <- xgb.importance(model = fit_best, feature_names = colnames(X))
cat("======== Top 20 most predictive features (by Gain) ========\n")
print(head(imp, 20))

# (Optional) To save importance as CSV:
# write.csv(imp, file = "xgb_feature_importance.csv", row.names = FALSE)

# (Optional) Quick plot of top-20:
# xgb.plot.importance(head(imp, 20))