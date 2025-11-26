library(BiocManager)
BiocManager::install("microbiome")
library(microbiome)
library(earth)
mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
source(paste0(mywd,'/Functions/HelperFunctions.R'))
data("dietswap")
ps <- dietswap

outcome <- estimate_richness(ps)$Shannon
otu_clr <- as.data.frame(t(otu_table(ps)))
fit <- earth(x = as.matrix(otu_clr), y = outcome)
y=outcome
X=otu_clr

pred <- predict(fit, newx = as.matrix(otu_clr), s="lambda.min")
cor(outcome, pred)^2  # R², usually >0.75
n <- nrow(otu_clr); p <- ncol(otu_clr)
data <- data.frame(cbind(y, otu_clr))
colnames(data) <- c("y", paste0("X", 1:p))
amountTrain <- 0.5
amountTest  <- 1 - amountTrain

# ---- parallel backend ----
n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1)
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
permR2Mars<-function(data,Y,j,model){
  Xperm <- data
  Xperm[, j] <- sample(data[, j], replace = FALSE)
  pred_perm <- predict(model, newdata = as.data.frame(Xperm))
  rsq_perm <- 1 - sum((Y - pred_perm)^2) / sum((Y - mean(Y))^2)
  rsq_perm
  return(rsq_perm)
}
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
                     dataTrain <- data[train_index, , drop = FALSE]
                     
                     mars_poly= earth(
                       y ~ .,
                       pmethod = "seqrep",
                       thresh=0,
                       data    = dataTrain
                     )
                     lm <- mars_poly
                     lm
                     # --- R^2 on train (not stored) / test halves (stored) ---?earth
                     pred1 <- predict(lm, newdata = as.data.frame(X[sample_index1, , drop = FALSE]))
                     pred2 <- predict(lm, newdata = as.data.frame(X[sample_index2, , drop = FALSE]))
                     y1 <- y[sample_index1]; y2 <- y[sample_index2]
                     
                     R2orig1 <- 1 - sum((y1 - pred1)^2) / sum((y1 - mean(y1))^2)
                     R2orig2 <- 1 - sum((y2 - pred2)^2) / sum((y2 - mean(y2))^2)
                     R2orig1;R2orig2
                     # --- permutation-based drops ---
                     Rnew1 <- sapply(seq_len(p), function(j)
                       permR2Mars(as.data.frame(X[sample_index1, , drop = FALSE]), Y = y1, j = j, model = lm))
                     Rnew2 <- sapply(seq_len(p), function(j)
                       permR2Mars(as.data.frame(X[sample_index2, , drop = FALSE]), Y = y2, j = j, model = lm))
                     
                     beta1 <- R2orig1 - Rnew1
                     beta2 <- R2orig2 - Rnew2
                     mirror <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
                     #hist(mirror)
                     selected_index <- SelectFeatures(mirror, abs(mirror),q=0.1)
                     num_sel <- length(selected_index)
                     num_sel
                     inc_row <- numeric(p)
                     fdp_val <- 0
                     pow_val <- 0
                     
                     if (num_sel > 0) {
                       inc_row[selected_index] <- 1 / num_sel
                     }
                     c(num_sel, R2orig1, R2orig2, inc_row)
                   }

parallel::stopCluster(cl)

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
}
selected_index
mean(R2orig1_vec)
mean(R2orig2_vec)