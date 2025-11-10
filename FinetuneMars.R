#Finetune Mars for the high-dimensional case



finetune=function(myseed,mynk, myfastk, myfastbeta, myminspan){
set.seed(myseed)
n <-400
p <- 500
p0 <- 25
q <- 0.1
delta <- 10
signal_index <- sample(c(1:p), size = p0, replace = F)
# simulate data
n1 <- floor(n/2); n2 <- n - n1
X1 <- matrix(rnorm(n1*p, mean= 0), n1, p)
X2 <- matrix(rnorm(n2*p, mean= 0), n2, p)
X  <- rbind(X1, X2)
beta_star <- numeric(p)
beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))*100
y <- (X^2 %*% beta_star + rnorm(n))
amountTrain <- 0.5
amountTest  <- 1 - amountTrain

n <- nrow(X); p <- ncol(X)
data <- data.frame(cbind(y, X))
colnames(data) <- c("y", paste0("X", 1:p))
size_half <- floor((amountTest/2) * n)
train_index <- sample.int(n, size = floor(amountTrain * n), replace = FALSE)
remaining_index <- setdiff(seq_len(n), train_index)
sample_index1 <- sample(remaining_index, size = size_half, replace = FALSE)
sample_index2 <- setdiff(remaining_index, sample_index1)
dataTrain <- data[train_index, , drop = FALSE]
# --- fit MARS ---
mars_poly <- earth(
  y ~ .,
  data    = dataTrain,
  degree  = 2,
  nk      = mynk,
  pmethod = "cv",
  nfold   = 5,
  ncross  = 3,
  trace   = 0,
  fast.k=myfastk,
  fast.beta=myfastbeta,
  minspan=myminspan
)
lm <- mars_poly
pred1 <- predict(lm, newdata = as.data.frame(X[sample_index1, , drop = FALSE]))
pred2 <- predict(lm, newdata = as.data.frame(X[sample_index2, , drop = FALSE]))
y1 <- y[sample_index1]; y2 <- y[sample_index2]

R2orig1 <- 1 - sum((y1 - pred1)^2) / sum((y1 - mean(y1))^2)
R2orig2 <- 1 - sum((y2 - pred2)^2) / sum((y2 - mean(y2))^2)

result=mean(R2orig1,R2orig2)
return(result)
}
# Define grid of (nk, seed)
nk_grid <- 30  # prepicked
n_seeds <- 3
myfast.k_grid=c(0,5,10,20)
myfast.beta_grid=c(0,1)
myminspan_grid=c(-3,0,1)
grid <- expand.grid(nk = nk_grid, seed = seq_len(n_seeds),myfastk=myfast.k_grid, myfastbeta=myfast.beta_grid, myminspan=myminspan_grid)

# Parallel setup
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
registerDoRNG(12345)

## Parallel grid search: for each nk, run 10 seeds i=1
res <- foreach(i = 1:nrow(grid),
               .combine = rbind,
               .inorder = FALSE,
               .packages = c("earth")) %dopar% {
                 nk <- grid$nk[i]
                 s  <- grid$seed[i]
                 myfastks  <- grid$myfastk[i]
                 myfastbetas  <- grid$myfastbeta[i]
                 myminspans  <- grid$myminspan[i]
                 r2 <- finetune(myseed = s, mynk = nk,
                                myfastk=myfastks, myfastbeta=myfastbetas, myminspan=myminspans)
                 data.frame(nk = nk, seed = s, myfastk=myfastks, myfastbeta=myfastbetas, myminspan=myminspans, R2 = r2)
               }

stopCluster(cl)

parallel::stopCluster(cl)

## Summarise across seeds for each nk
summary_res <- res %>%
  group_by(myfastk, myfastbeta, myminspan) %>%
  summarise(
    mean_R2 = mean(R2),
    sd_R2   = sd(R2),
    q25     = quantile(R2, 0.25),
    q75     = quantile(R2, 0.75),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_R2))

print(summary_res, n = nrow(summary_res))

best_nk    <- summary_res$nk[1]
best_stats <- summary_res[1, ]

cat(sprintf(
  "Best nk = %d | mean R^2 = %.4f (sd = %.4f; IQR = [%.4f, %.4f])\n",
  best_nk, best_stats$mean_R2, best_stats$sd_R2, best_stats$q25, best_stats$q75
))

## (Optional) Fit once at best nk with a fixed seed to get a final score:
## set.seed(999)
## final_score <- finetune(myseed = 999, mynk = best_nk)
## cat(sprintf("Final check at nk=%d with seed=999: R^2 = %.4f\n", best_nk, final_score))