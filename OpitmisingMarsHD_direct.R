#Finetune Mars for the high-dimensional case

library(earth)
source('C:/Users/mde4023/Downloads/FDR_Datasplitting/Functions/MarsParallelHD.R')
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
    myprmethod=myprmethod
) {
  # Derive default per-stage seeds (stable & independent) if not provided
  # --- choose signal indices ---
  set.seed(myseed)
  signal_index <- sample.int(p, size = p0, replace = FALSE)
  
  # --- simulate data ---
  n1 <- floor(n/2); n2 <- n - n1
  X1 <- matrix(rnorm(n1 * p, mean = -1), n1, p)
  X2 <- matrix(rnorm(n2 * p, mean = 1), n2, p)
  X  <- rbind(X1, X2)
  
  beta_star <- numeric(p)
  beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta * sqrt(log(p) / n)) * 100
  
  # Nonlinear quadratic signal + Gaussian noise
  y <- as.vector((X^2 %*% beta_star) + rnorm(n, sd = noise_sd))
  
  # --- fit MARS --- ?earth
  g1 <- ApplyMarsTrain_HDparallel( X = X, y = y, q = 0.10, num_split = 5,mynk=mynk, ,mypmethod=mypmethod,myprmethod=myprmethod, signal_index = signal_index, myseed = 1)
  
  return(g1$MDS_power)
}

# Define grid of (nk, seed)

nk_grid <- c(30,40,50,60,70,80,90,100,110,120) # adjust range
mypmethod_grid   <- c("backward", "seqrep", "cv", "forward")
myprmethod_g=c(5,10,25)

grid_full <- expand.grid(
  seed=c(11132025,11132026,11132027),
  nk        = nk_grid,
  pmethod   = mypmethod_grid,
  myprmethod=myprmethod_g
)
result=c()
for(i in 1:nrow(grid_full)){
  seed  <- grid_full$seed[i]
  nk      <- grid_full$nk[i]
  pmethod <- as.character(grid_full$pmethod[i])
  myprmethod=grid_full$myprmethod[i]
  
  r2 <- finetune(
    myseed      = seed,
    mynk        = nk,
    mypmethod   = pmethod,
    myprmethod=myprmethod
  )
  nresult=data.frame(
    seed=seed,
    nk        = nk,
    pmethod   = pmethod,
    myprmethod=myprmethod,
    MDS_power        = r2
  )
  result=rbind(result,nresult)
  print(result)
}
res <- foreach(i = seq_len(nrow(grid_full)),
               .combine = rbind,
               .inorder = FALSE,
               .packages = c("earth")) %dopar% {
                 seed  <- grid_full$seed[i]
                 nk      <- grid_full$nk[i]
                 pmethod <- grid_full$pmethod[i]
                 nprune  <- grid_full$nprune[i]
                 thresh  <- grid_full$thresh[i]
                 penalty <- grid_full$penalty[i]
                 fast.k  <- grid_full$fast.k[i]
                 fast.beta <- grid_full$fast.beta[i]
                 minspan <- grid_full$minspan[i]
                 
                 r2 <- finetune(
                   myseed      = seed,
                   mynk        = nk,
                   mypmethod   = pmethod,
                   mynprune    = nprune,
                   mypenalty   = penalty,
                   mythresh =thresh,
                   myfast.k    = fast.k,
                   myfast.beta = fast.beta,
                   myminspan   = minspan
                 )
                 
                 # INTERMEDIATE OUTPUT (appears in console)
                 cat(sprintf(
                   "Done %d/%d | nk=%d | pmethod=%s | R2=%.4f\n",
                   i, nrow(grid_full), nk, pmethod, r2
                 ))
                 
                 data.frame(
                   nk        = nk,
                   pmethod   = pmethod,
                   nprune    = nprune,
                     penalty   = penalty,
                   fast.k    = fast.k,
                   thresh =thresh,
                   fast.beta = fast.beta,
                   minspan   = minspan,
                   R2        = r2
                 )
               }
stopCluster(cl)
View(res)
summary_res <- res %>%
  group_by(nk, pmethod, nprune,thresh, penalty, fast.k, fast.beta, minspan) %>%
  summarise(
    mean_R2 = mean(R2),
    sd_R2   = sd(R2),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_R2))
newgrid=summary_res[which(round(summary_res$mean_R2,3)==round(max(summary_res$mean_R2),3)),]

seeds <- 121:150   # the 10 seeds you want to use

# newgrid has 12 rows in your example
n_rows <- nrow(newgrid)

# Repeat each row 10 times and add a seed column
grid_full <- newgrid[rep(seq_len(n_rows), each = length(seeds)), ]
grid_full$seed <- rep(seeds, times = n_rows)

# (Optional) drop old R2 column if you don't want it in the run
grid_full$R2 <- NULL

head(grid_full)
nrow(grid_full) 
suppressPackageStartupMessages({
  library(earth)
  library(foreach)
  library(doParallel)
  library(doRNG)
  library(dplyr)
})

# Parallel setup
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
registerDoRNG(12345)


res <- foreach(i = seq_len(nrow(grid_full)),
               .combine = rbind,
               .inorder = FALSE,
               .packages = "earth") %dopar% {
                 
                 nk        <- grid_full$nk[i]
                 pmethod   <- grid_full$pmethod[i]
                 nprune    <- grid_full$nprune[i]
                 penalty   <- grid_full$penalty[i]
                 fast.k    <- grid_full$fast.k[i]
                 thresh  <- grid_full$thresh[i]
                 fast.beta <- grid_full$fast.beta[i]
                 minspan   <- grid_full$minspan[i]
                 seed      <- grid_full$seed[i]
                 
                 r2 <- finetune(
                   myseed      = seed,
                   mynk        = nk,
                   mypmethod   = pmethod,
                   mynprune    = nprune,
                   mypenalty   = penalty,
                   mythresh =thresh,
                   myfast.k    = fast.k,
                   myfast.beta = fast.beta,
                   myminspan   = minspan
                 )
                 
                 # ---- INTERMEDIATE OUTPUT ----
                 cat(sprintf(
                   "Done %3d/%3d | seed=%2d | nk=%d | pmethod=%s | penalty=%d | fast.k=%d | fast.beta=%d | minspan=%d | R2=%.4f\n",
                   i, nrow(grid_full), seed, nk, pmethod, penalty, fast.k, fast.beta, minspan, r2
                 ))
                 
                 data.frame(
                   nk        = nk,
                   pmethod   = pmethod,
                   nprune    = nprune,
                   penalty   = penalty,
                   fast.k    = fast.k,
                   fast.beta = fast.beta,
                   minspan   = minspan,
                   thresh =thresh,
                   seed      = seed,
                   R2        = r2
                 )
               }

stopCluster(cl)      

summary_res <- res %>%
  group_by(nk, pmethod, nprune, penalty, fast.k, fast.beta, minspan) %>%
  summarise(
    mean_R2 = mean(R2),
    sd_R2   = sd(R2),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_R2))

print(summary_res, n = nrow(summary_res))
