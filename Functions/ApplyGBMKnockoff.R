ApplyGBMKnockoff <- function(X, y, q = 0.1,
                             nrounds = 500,
                             param =list(
                               objective = "reg:squarederror",
                               eta       = 0.05,
                               max_depth = 3,
                               subsample = 0.6,
                               colsample_bytree = 0.8,
                               lambda    = 1,
                               alpha     = 0
                             ),
                             seed = 1) {
  
  set.seed(seed)
  
  ## 1. build knockoff copies (second-order Gaussian)
  ko <- create.second_order(X)             # returns list with X_k, Sigma, etc.
  Xk <- ko
  
  ## 2. fit XGBoost on [X, X_k]
  X_aug <- cbind(X, Xk)
  colnames(X_aug) <- c(paste0("X",  seq_len(ncol(X))),
                       paste0("Xk", seq_len(ncol(X))))
  dtrain <- xgb.DMatrix(data = X_aug, label = y)
  
  bst <- xgb.train(params   = param,
                   data     = dtrain,
                   nrounds  = nrounds,
                   verbose  = 0)
  
  ## 3. feature-gain importance for originals & knockoffs
  imp <- xgb.importance(model = bst, feature_names = colnames(X_aug))
  gain <- numeric(ncol(X_aug)); names(gain) <- colnames(X_aug)
  gain[imp$Feature] <- imp$Gain            # features not used keep gain = 0
  
  ## 4. antisymmetric statistic
  p <- ncol(X)
  W <- gain[1:p] - gain[(p + 1):(2 * p)]
  names(W) <- paste0("X", seq_len(p))
  
  ## 5. apply knockoff threshold
  T <- knockoff.threshold(W, fdr = q, offset = 1)   # BC+1 threshold
  selected_index <- which(W >= T)
  result <- CalculateFDP_Power(selected_index, signal_index)
  MDS_fdp <- result$fdp
  MDS_power <- result$power
  list(fdp  = MDS_fdp,           # indices of discoveries
       power         = MDS_power)
}
