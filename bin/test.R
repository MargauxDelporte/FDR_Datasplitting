res_mat <- foreach(iter = 1:2,
                   .combine = "rbind",
                   .packages = c("randomForest")) %dorng% {
                     
                     permR2 <- function(data, Y, j, model) {
                       Xperm <- data
                       # Permute column j
                       Xperm[, j] <- sample(data[, j], replace = FALSE)
                       
                       # Predict using permuted data
                       pred_perm <- predict(model, newdata = as.data.frame(Xperm))
                       
                       # Calculate R2
                       rsq_perm <- 1 - sum((Y - pred_perm)^2) / sum((Y - mean(Y))^2)
                       return(rsq_perm)
                     }
                     source(paste0('C:/Users/mde4023/Downloads/FDR_Datasplitting','/Functions/HelperFunctions.R'))
                     p=1043;n=130
                     data_full=mydata_full
                     X=mydata_full[,-1]
                     y=mydata_full[,1]
                     # --- indices ---
                     train_index     <- sample.int(n, size = floor(amountTrain * n), replace = FALSE)
                     remaining_index <- setdiff(seq_len(n), train_index)
                     
                     # split the remaining part in two halves
                     size_half     <- floor((amountTest / 2) * n)
                     sample_index1 <- sample(remaining_index, size = size_half, replace = FALSE)
                     sample_index2 <- setdiff(remaining_index, sample_index1)
                     
                     dataTrain <- data_full[train_index, , drop = FALSE]
                     
                     # --- fit RF using parameter vector pm ---
                     mynlm <- randomForest(
                       y ~ ., 
                       ntry=100,
                       ntree=104,
                       nodesize=1,
                       data    = data_full
                     )
                     
                     # --- R² on the two halves ---
                     pred1 <- predict(mynlm, newdata = as.data.frame(X[sample_index1, , drop = FALSE]))
                     pred2 <- predict(mynlm, newdata = as.data.frame(X[sample_index2, , drop = FALSE]))
                     
                     y1 <- y[sample_index1]
                     y2 <- y[sample_index2]
                     
                     R2orig1 <- 1 - sum((y1 - pred1)^2) / sum((y1 - mean(y1))^2)
                     R2orig2 <- 1 - sum((y2 - pred2)^2) / sum((y2 - mean(y2))^2)
                     
                     # --- permutation-based drops ---
                     Rnew1 <- sapply(seq_len(p), function(j)
                       permR2(as.data.frame(X[sample_index1, , drop = FALSE]), Y = y1, j = j, model = mynlm))
                     Rnew2 <- sapply(seq_len(p), function(j)
                       permR2(as.data.frame(X[sample_index2, , drop = FALSE]), Y = y2, j = j, model = mynlm))
                     
                     beta1  <- R2orig1 - Rnew1
                     beta2  <- R2orig2 - Rnew2
                     mirror <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
                     
                     selected_index <- SelectFeatures(mirror, abs(mirror), q = 0.1)
                     num_sel <- length(selected_index)
                     
                     inc_row <- numeric(p)
                     if (num_sel > 0) {
                       inc_row[selected_index] <- 1 / num_sel
                     }
                     
                     c(num_sel, R2orig1, R2orig2, inc_row)
                   }

# ---- unpack ----
num_select     <- res_mat[, 1]
R2orig1_vec    <- res_mat[, 2]
R2orig2_vec    <- res_mat[, 3]
inclusion_rate_mat <- res_mat[, -(1:3), drop = FALSE]  # num_split x p

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
}else(
  selected_index=c()
)
selected_index
