# ------------------------- Setup -------------------------
library(softbart)   # install.packages("softbart") if needed

set.seed(42)

R2 <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)

# Quadratic screening (train-only)
library(softbart)
# helper: pass only args supported by your installed softbart::softbart()
sb_fun <- SoftBart::softbart
sb_formals <- names(formals(sb_fun))
mk_args <- function(common_args, cand_args) {
  args <- c(common_args, cand_args)
  args[names(args) %in% sb_formals]
}

# extract posterior-mean predictions from a softbart fit
pred_mean_softbart <- function(fit, Xtest) {
  if (!is.null(fit$yhat.test.mean)) return(as.vector(fit$yhat.test.mean))
  if (!is.null(fit$yhat.test))      return(colMeans(fit$yhat.test))
  # fall back to predict() if available
  if ("predict" %in% methods(class = class(fit))) {
    p <- predict(fit, newdata = Xtest)
    if (is.matrix(p)) return(colMeans(p))
    return(as.vector(p))
  }
  stop("Could not extract predictions from softbart fit.")
}

# ------------------------- Data prep -------------------------
X_train <- as.matrix(X_train); X_test <- as.matrix(X_test)
stopifnot(all(is.finite(X_train)), all(is.finite(X_test)))
# drop any zero-variance cols (train-defined)
keep_cols <- 1:p
X_train <- X_train[, keep_cols, drop = FALSE]
X_test  <- X_test[,  keep_cols, drop = FALSE]
Xtr=X_train
Xte=X_test
# scale y for stability
y_mu <- mean(y_train); y_sd <- sd(y_train); stopifnot(y_sd > 0)
ytr  <- (y_train - y_mu) / y_sd
y_t<- scale(y_test)
# ------------------------- CV grid -------------------------
# Notes:
# - alpha/beta: tree-depth prior (larger beta => shallower trees)
# - tau: softness of splits (larger => more "hard" BART-like)
# - sparse/theta: Dirichlet prior for variable selection; theta ≈ expected # active vars
grid <- expand.grid(
  ntree  = c(50, 100, 200),
  alpha  = c(0.95, 0.85),
  beta   = c(2, 1),
  tau    = c(1, 3, 6),
  sparse = TRUE,
  theta  = c(10, 25),   # p0 ~ 25 in your setup; tune around it
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)

K <- 5  # 5-fold CV on the (screened) training set
fold_id <- sample(rep(1:K, length.out = nrow(Xtr)))

# sampling params (keep modest for speed; increase for final runs)
nskip  <- 800
ndpost <- 800

cv_results <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  r2_fold <- numeric(K)
  
  for (k in 1:K) {
    tr <- which(fold_id != k)
    va <- which(fold_id == k)
    
    common <- list(
      x.train = Xtr[tr, , drop = FALSE],
      y.train = ytr[tr],
      x.test  = Xtr[va, , drop = FALSE],
      nskip   = nskip,
      ndpost  = ndpost
    )
    hypers <- Hypers(X = X_train, Y = ytr,
                     num_tree = 50, beta = 1, gamma = 0.9)
    opts <- Opts(num_burn = 5000, num_save = 5000,
                 update_s = FALSE, update_sigma_mu = FALSE)
    
    # build args that your version supports
    fitted_softbart_2 <- softbart(X_train, ytr, Xte,
                                  opts = opts, hypers = hypers)
    R2(fitted_softbart_2$y_hat_test_mean, y_t)
    variable_selection <- posterior_probs(fitted_softbart_2)
    plot(variable_selection$post_probs,
         col = ifelse(1:250 < 6, "#386CB0", "#7FC97F"), pch = 20,
         xlab = "Predictor", ylab = "PIP", main = "Variable Selection",
         ylim = c(0,1))
    abline(h = 0.5, col = "darkblue", lwd = 2, lty = 3)
    which(variable_selection$post_probs>0.25)%in%signal_index
    if (inherits(fit, "try-error")) { r2_fold[k] <- NA; next }
    
    yhat_val_scaled <- pred_mean_softbart(fit, Xtr[va, , drop = FALSE])
    yhat_val <- y_mu + y_sd * yhat_val_scaled
    r2_fold[k] <- R2(y_train[va], yhat_val)
  }
  
  cv_results[[i]] <- data.frame(
    i = i,
    mean_R2 = mean(r2_fold, na.rm = TRUE),
    sd_R2   = sd(r2_fold,  na.rm = TRUE),
    fail    = sum(is.na(r2_fold))
  )
}

cv_tab <- do.call(rbind, cv_results)
best_i <- with(cv_tab, i[which.max(mean_R2)])
best_hypers <- grid[best_i, , drop = FALSE]

cat("Best (CV):\n"); print(cbind(best_hypers, cv_tab[cv_tab$i == best_i, -1]))

# ------------------------- Refit on all training & test -------------------------
common_full <- list(
  x.train = Xtr, y.train = ytr,
  x.test  = Xte,
  nskip   = max(nskip, 1500),   # beef up for final fit
  ndpost  = max(ndpost, 1500)
)
fit_best <- do.call(sb_fun, mk_args(common_full, as.list(best_hypers)))

pred_test_scaled <- pred_mean_softbart(fit_best, Xte)
pred_test <- y_mu + y_sd * pred_test_scaled
R2_test <- R2(y_test, pred_test)
cat(sprintf("SoftBART (screened m=%d) | Test R^2 = %.3f\n", length(idx), R2_test))