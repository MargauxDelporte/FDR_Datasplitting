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


result[which(result$MDS_power>0.5),]
#> result[which(result$MDS_power>0.5),]
#seed  nk  pmethod myprmethod MDS_power
#245 11132026  40 backward         25      0.52
#248 11132026  50 backward         25      0.52
#251 11132026  60 backward         25      0.52
#254 11132026  70 backward         25      0.52
#257 11132026  80 backward         25      0.52
#260 11132026  90 backward         25      0.52
#263 11132026 100 backward         25      0.52
#266 11132026 110 backward         25      0.52
#269 11132026 120 backward         25      0.52

nk_grid <- seq(from=30,to=150,by=10) # adjust range
mypmethod_grid   <- c("backward")
myprmethod_g=seq(from=25,to=35,by=5)

grid_full <- expand.grid(
  seed=c(11142025,11142026,11142027),
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
result[which(result$MDS_power>0.47),]

#> result[which(result$MDS_power>0.47),]
#seed  nk  pmethod myprmethod MDS_power
#46 11142025  50 backward         30      0.48
#49 11142025  60 backward         30      0.48
#52 11142025  70 backward         30      0.48
#55 11142025  80 backward         30      0.48
#58 11142025  90 backward         30      0.48
#61 11142025 100 backward         30      0.48
#64 11142025 110 backward         30      0.48
#67 11142025 120 backward         30      0.48
#70 11142025 130 backward         30      0.48
#73 11142025 140 backward         30      0.48
#76 11142025 150 backward         30      0.48