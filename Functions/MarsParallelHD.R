suppressPackageStartupMessages({
  library(earth)        # MARS
  library(foreach)
  library(doParallel)
  library(doRNG)        # reproducible %dopar%
})


#setwd("~/Desktop/temp") j=1 model=lm

permR2Mars<-function(data,Y,j,model){
  Xperm <- data
  Xperm[, j] <- sample(data[, j], replace = FALSE)
  pred_perm <- predict(model, newdata = as.data.frame(Xperm))
  rsq_perm <- 1 - sum((Y - pred_perm)^2) / sum((Y - mean(Y))^2)
  rsq_perm
  return(rsq_perm)
}
ApplyMarsTrain_parallel <- function(X, y, q, myseed, num_split = 50,
                                    signal_index = signal_index,
                                    plot_hist = FALSE) {
  
  stopifnot(nrow(X) == length(y))
  set.seed(myseed)
  
  amountTrain <- 0.5
  amountTest  <- 1 - amountTrain
  
  n <- nrow(X); p <- ncol(X)
  data <- data.frame(cbind(y, X))
  colnames(data) <- c("y", paste0("X", 1:p))
  
  # ---- parallel backend ----
  n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1)
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  # We’ll collect per-split:
  # [1] num_selected, [2] fdp, [3] power,
  # [4] R2orig1, [5] R2orig2, [6:(5+p)] inclusion_rate_row (length p)
  res_mat <- foreach(iter = 1:num_split,
                     .combine = "rbind",
                     .packages = c("earth"),
                     .export   = c("permR2Mars","SelectFeatures","CalculateFDP_Power")) %dorng% {
                       
                       # --- indices ---
                       train_index <- sample.int(n, size = floor(amountTrain * n), replace = FALSE)
                       remaining_index <- setdiff(seq_len(n), train_index)
                       
                       # split the remaining half evenly
                       size_half <- floor((amountTest/2) * n)
                       sample_index1 <- sample(remaining_index, size = size_half, replace = FALSE)
                       sample_index2 <- setdiff(remaining_index, sample_index1)
                       
                       # --- fit MARS ---
                       dataTrain <- data[train_index, , drop = FALSE]
                       mars_poly <- earth(
                         y ~ .,
                         data    = dataTrain,
                         degree  = 2,
                         nk      = 30,
                         pmethod = "cv",
                         nfold   = 5,
                         ncross  = 3,
                         trace   = 0,
                         fast.k=5,
                         fast.beta=1,
                         minspan=-3
                       )
                       lm <- mars_poly
                       
                       # --- R^2 on train (not stored) / test halves (stored) ---
                       pred1 <- predict(lm, newdata = as.data.frame(X[sample_index1, , drop = FALSE]))
                       pred2 <- predict(lm, newdata = as.data.frame(X[sample_index2, , drop = FALSE]))
                       y1 <- y[sample_index1]; y2 <- y[sample_index2]
                       
                       R2orig1 <- 1 - sum((y1 - pred1)^2) / sum((y1 - mean(y1))^2)
                       R2orig2 <- 1 - sum((y2 - pred2)^2) / sum((y2 - mean(y2))^2)
                       
                       # --- permutation-based drops ---
                       Rnew1 <- sapply(seq_len(p), function(j)
                         permR2Mars(as.data.frame(X[sample_index1, , drop = FALSE]), Y = y1, j = j, model = lm))
                       Rnew2 <- sapply(seq_len(p), function(j)
                         permR2Mars(as.data.frame(X[sample_index2, , drop = FALSE]), Y = y2, j = j, model = lm))
                       
                       beta1 <- R2orig1 - Rnew1
                       beta2 <- R2orig2 - Rnew2
                       mirror <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
                       
                       # optional: plotting is off by default in parallel
                       if (plot_hist && length(setdiff(seq_len(p), signal_index)) > 1) {
                         # store stats silently, or skip to avoid device contention
                         # hist(mirror[-signal_index])  # not recommended inside parallel loop
                         invisible(NULL)
                       }
                       
                       selected_index <- SelectFeatures(mirror, abs(mirror), q)
                       
                       num_sel <- length(selected_index)
                       inc_row <- numeric(p)
                       fdp_val <- NA_real_
                       pow_val <- NA_real_
                       
                       if (num_sel > 0) {
                         inc_row[selected_index] <- 1 / num_sel
                         res <- CalculateFDP_Power(selected_index, signal_index)
                         fdp_val <- res$fdp
                         pow_val <- res$power
                       }
                       
                       c(num_sel, fdp_val, pow_val, R2orig1, R2orig2, inc_row)
                     }
  
  parallel::stopCluster(cl)
  
  # ---- unpack ----
  num_select     <- res_mat[, 1]
  fdp            <- res_mat[, 2]
  power          <- res_mat[, 3]
  R2orig1_vec    <- res_mat[, 4]
  R2orig2_vec    <- res_mat[, 5]
  inclusion_rate_mat <- res_mat[, -(1:5), drop = FALSE]  # num_split x p
  
  # DS results from the *first* split (order preserved with foreach)
  DS_fdp   <- fdp[1]
  DS_power <- power[1]
  
  # MDS inclusion rates (avg across splits)
  inclusion_rate <- apply(inclusion_rate_mat, 2, mean)
  
  # Rank features by empirical inclusion rate
  feature_rank <- order(inclusion_rate)
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  
  if (length(feature_rank) != 0) {
    null_feature <- numeric()
    for (feature_index in seq_along(feature_rank)) {
      if (sum(inclusion_rate[feature_rank[1:feature_index]]) > q) break
      null_feature <- c(null_feature, feature_rank[feature_index])
    }
    selected_index <- setdiff(feature_rank, null_feature)
    
    res_final <- CalculateFDP_Power(selected_index, signal_index)
    MDS_fdp   <- res_final$fdp
    MDS_power <- res_final$power
  } else {
    MDS_fdp <- 0
    MDS_power <- 0
  }
  
  message(paste0("First R squared: ", round(R2orig1_vec[1], 3)))
    message(paste0("Second R squared: ", round(R2orig2_vec[1], 3)))
    message(paste0("DS_fdp = ", DS_fdp,
                " DS_power = ", DS_power,
                                 " MDS_fdp = ", MDS_fdp,
                                 " MDS_power = ", MDS_power))
  
  return(list(
    DS_fdp   = ifelse(is.na(DS_fdp),0,DS_fdp),
    DS_power = ifelse(is.na(DS_power),0,DS_power),
    MDS_fdp  = MDS_fdp,
    MDS_power = MDS_power)
  )
}



