#Finetune Mars for the high-dimensional case



finetune <- function(
    myseed, mynk,
    # --- data generation ---
    n = 400, p = 500, p0 = 25,              # sample size, total predictors, #signals
    delta = 10,                             # signal scale (on beta)
    noise_sd = 1,                           # sd of N(0, noise_sd^2) noise
    # --- split ---
    amountTrain = 0.5,                      # training fraction
    # --- MARS (earth) options ---
    penalty = 2,
    mypmethod = "cv",
    nfold = 5,
    myfast.k = 50,
    trace = 0,
    myfast.beta=1,
    myminspan,
    mythresh,
    mynprune,
    mypenalty
) {
  # Derive default per-stage seeds (stable & independent) if not provided
  # --- choose signal indices ---
  set.seed(myseed)
  signal_index <- sample.int(p, size = p0, replace = FALSE)
  
  # --- simulate data ---
  n1 <- floor(n/2); n2 <- n - n1
  set.seed(data_seed)
  X1 <- matrix(rnorm(n1 * p, mean = 0), n1, p)
  X2 <- matrix(rnorm(n2 * p, mean = 0), n2, p)
  X  <- rbind(X1, X2)
  
  beta_star <- numeric(p)
  beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta * sqrt(log(p) / n)) * 100
  
  # Nonlinear quadratic signal + Gaussian noise
  y <- as.vector((X^2 %*% beta_star) + rnorm(n, sd = noise_sd))
  
  # --- prepare frames & split ---
  amountTest <- 1 - amountTrain
  data <- data.frame(y = y, X)
  names(data) <- c("y", paste0("X", seq_len(p)))
  
  size_half <- floor((amountTest / 2) * n)
  
  set.seed(split_seed)
  train_index     <- sample.int(n, size = floor(amountTrain * n), replace = FALSE)
  remaining_index <- setdiff(seq_len(n), train_index)
  sample_index1   <- sample(remaining_index, size = size_half, replace = FALSE)
  sample_index2   <- setdiff(remaining_index, sample_index1)
  
  dataTrain <- data[train_index, , drop = FALSE]
  
  # --- fit MARS ---
  mars_fit <- earth(
    y ~ .,
    data    = dataTrain,
    degree  = 2,
    nk      = mynk,
    pmethod = mypmethod,
    nfold   = 5,
    ncross  = 3,
    penalty   = mypenalty,
    fast.k= myfast.k,
    fast.beta= myfast.beta,
    minspan=myminspan,
   thresh=mythresh,
   nprune=mynprune
  )
  
  # --- evaluate on two disjoint test halves ---
  pred1 <- predict(mars_fit, newdata = data[sample_index1, , drop = FALSE])
  pred2 <- predict(mars_fit, newdata = data[sample_index2, , drop = FALSE])
  
  y1 <- y[sample_index1]; y2 <- y[sample_index2]
  
  R2orig1 <- 1 - sum((y1 - pred1)^2) / sum((y1 - mean(y1))^2)
  R2orig2 <- 1 - sum((y2 - pred2)^2) / sum((y2 - mean(y2))^2)
  
  mean(c(R2orig1, R2orig2))
}

# Define grid of (nk, seed)
nk_grid <- seq(10, 50, by = 10)  # adjust range
mypmethod_grid=c('backward', 'none', 'exhaustive', 'forward', 'seqrep', 'cv')
mynprune_grid=c(10,25,40)
mythresh_grid=c(0.0001,0.001,0.01)
mypenalty_grid=c(-1,0,2,3,4)
myfast.k_grid=c(0,5,10,20)
myfast.beta_grid=c(0,1)
myminspan_grid=c(-3,0,1)

# Create full parameter grid
grid <- expand.grid(
  nk         = nk_grid,
  pmethod    = mypmethod_grid,
  nprune     = mynprune_grid,
  thresh     = mythresh_grid,
  penalty    = mypenalty_grid,
  fast.k     = myfast.k_grid,
  fast.beta  = myfast.beta_grid,
  minspan    = myminspan_grid,
  myseed       = 1:3,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)




















































suppressPackageStartupMessages({
  library(earth)
  library(dplyr)
  library(purrr)
})


# Preview grid size and first rows
nrow(grid)
head(grid)

# Parallel setup
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
registerDoRNG(12345)

## Parallel grid search: for each nk, run 10 seeds
res <- foreach(i = 1:nrow(grid),
               .combine = rbind,
               .inorder = FALSE,
               .packages = c("earth")) %dopar% {
                 nk <- grid$nk[i]
                 s  <- grid$seed[i]
                 r2 <- finetune(myseed = s, mynk = nk)
                 data.frame(nk = nk, seed = s, R2 = r2)
               }

stopCluster(cl)
warnings() 
parallel::stopCluster(cl)
library(dplyr)
## Summarise across seeds for each nk
summary_res <- res %>%
  group_by(nk) %>%
  summarise(
    mean_R2 = mean(R2),
    sd_R2   = sd(R2),
    q25     = quantile(R2, 0.25),
    q75     = quantile(R2, 0.75),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_R2))

print(summary_res, n = nrow(summary_res))
View(res)
best_nk    <- summary_res$nk[1]
best_stats <- summary_res[1, ]

cat(sprintf(
  "Best nk = %d | mean R^2 = %.4f (sd = %.4f; IQR = [%.4f, %.4f])\n",
  best_nk, best_stats$mean_R2, best_stats$sd_R2, best_stats$q25, best_stats$q75
))

finetune(8,30)
?earth
## (Optional) Fit once at best nk with a fixed seed to get a final score:
## set.seed(999)
## final_score <- finetune(myseed = 999, mynk = best_nk)
## cat(sprintf("Final check at nk=%d with seed=999: R^2 = %.4f\n", best_nk, final_score))