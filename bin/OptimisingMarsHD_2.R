## Install needed packages if you don't have them
## install.packages("earth")
## install.packages("glmnet")

library(earth)   # MARS
library(glmnet)  # ridge / lasso
source('C:/Users/mde4023/Downloads/FDR_Datasplitting/Functions/MarsParallelHD.R')
#----- 1. Single simulation run -------------------------------------------

sim_once <- function(n = 400, p = 500, sigma = 1) {
  set.seed(123)
  # Generate predictors: X_ij ~ Uniform(0,1)
  X <- matrix(runif(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", 1:p)
  
  # Random coefficients for the 25-signal nonlinear function
  a <- rnorm(5, mean = 3, sd = 1)
  b <- rnorm(5, mean = 2, sd = 1)
  c <- rnorm(5, mean = 3, sd = 1)
  d <- rnorm(5, mean = 2, sd = 1)
  
  f <- 0
  for (k in 0:4) {
    idx <- 5 * k
    f <- f +
      a[k + 1] * sin(pi * X[, idx + 1]) +
      b[k + 1] * (X[, idx + 2] - 0.5)^2 +
      c[k + 1] * X[, idx + 3] * X[, idx + 4] +
      d[k + 1] * pmax(0, X[, idx + 5] - 0.3)
  }
  # Add noise
  y <- 1000*f + rnorm(n, mean = 0, sd = sigma)
  signal_index=1:25
  g1 <- ApplyMarsTrain_HDparallel( X = X, y = y, q = 0.10, num_split = 5, signal_index = signal_index, myseed = 1)
  
  return(g1$MDS_power)
}
sim_once()

sim_once <- function(mnk,mprune,seed) {
  n=400
  p=500
  set.seed(seed)
  p0=25
  delta=10
  noise_sd=1
  signal_index <- sample.int(p, size = 25, replace = FALSE)
  
  # --- simulate data ---
  n1 <- floor(n/2); n2 <- n - n1
  X1 <- matrix(rnorm(n1 * p, mean = -1), n1, p)
  X2 <- matrix(rnorm(n2 * p, mean = 1), n2, p)
  X  <- rbind(X1, X2)
  
  beta_star <- numeric(p)
  beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta * sqrt(log(p) / n)) * 1000
  
  # Nonlinear quadratic signal + Gaussian noise
  y <- as.vector((X^2 %*% beta_star) + rnorm(n, sd = noise_sd))
  g1 <- ApplyMarsTrain_HDparallel( X = X, y = y, q = 0.10, num_split = 5,mynk=mnk,myprune=mprune, signal_index = signal_index, myseed = 1)
  
  return(g1$MDS_power)
}

result=c()
for(i in seq(from=30,to=200,by=10)){
  for(j in seq(from=25,to=i,by=10)){
    for(k in c(11,123,569)){
    r2=sim_once(mnk=i,mprune=j,k)
    nresult=data.frame(
      seed=k,
      nk        = i,
      pruned   = j,
      MDS_power        = r2
    )
    result=rbind(result,nresult)
    print(result)
  }
  }
}
sim_once(mnk=50,mprune=30)
