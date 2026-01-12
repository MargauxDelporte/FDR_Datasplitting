### High dimension linear model
rm(list = ls())

mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
setwd(mywd)

source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Scenario/Scenario3b/TriangleLassoHD2.R'))
source(paste0(mywd,'/Scenario/Scenario3b/MarsParallelHD.R'))
source(paste0(mywd,'/Functions Dai/knockoff.R'))
source(paste0(mywd,'/Functions Dai/analysis.R'))
source(paste0(mywd,'/Functions Dai/MBHq.R'))
source(paste0(mywd,'/Functions Dai/DS.R'))
source(paste0(mywd,'/Functions Dai/fdp_power.R'))

#devtools::install_github("Jeremy690/DSfdr/DSfdr",force = TRUE)
library(MASS)
library(glmnet)
library(knockoff)
library(mvtnorm)
library(hdi)
library(parallel)

### algorithmic settings
num_split <- 5
n <-400
p <- 500
p0 <- 25
q <- 0.1
set.seed(456)
signal_index <- sample(c(1:p), size = p0, replace = F)

#######set up the method for the comparison############# 
S <- 5
seeds <- 1:S

lambda_grid <- seq(from = 0, to = 0.1, by = 0.0005)

# store everything (seed, lambda, fdp, power)
all_results <- data.frame(seed = integer(0),
                          lambda = double(0),
                          fdp = double(0),
                          power = double(0))

for (s in seeds) {
  set.seed(s)
  
  delta <- 10  # your signal strength
  
  # simulate data
  n1 <- floor(n/2); n2 <- n - n1
  X1 <- matrix(rnorm(n1*p, mean = -1), n1, p)
  X2 <- matrix(rnorm(n2*p, mean =  1), n2, p)
  X  <- rbind(X1, X2)
  
  beta_star <- numeric(p)
  beta_star[signal_index] <- rnorm(p0, 0, delta * sqrt(log(p)/n))
  y <- as.numeric(X %*% beta_star + rnorm(n))
  
  # run method over lambda grid
  for (lambda in lambda_grid) {
    g1 <- ApplyMarsTrain_HDparallel(
      X = X, y = y, q = q, num_split = 5,
      signal_index = signal_index,
      myseed = 1,  
    )
    
    all_results <- rbind(all_results,
                         data.frame(seed = s,
                                    lambda = lambda,
                                    fdp = g1$MDS_fdp,
                                    power = g1$MDS_power))
    print(lambda)
  }
}

# average over seeds for each lambda
avg_results <- aggregate(cbind(fdp, power) ~ lambda, data = all_results, FUN = mean)

# optional: also get SD over seeds
sd_results  <- aggregate(cbind(fdp, power) ~ lambda, data = all_results, FUN = sd)

head(avg_results)
 
View(avg_results)

myresults=c()
i=10
for(s in seq(from=1,to=10,by=01)){
  set.seed(s)
  myresults=c()
  delta <- i
  # simulate data
  n1 <- floor(n/2); n2 <- n - n1
  X1 <- matrix(rnorm(n1*p, mean= -1), n1, p)
  X2 <- matrix(rnorm(n2*p, mean= 1), n2, p)
  X  <- rbind(X1, X2)
  beta_star <- numeric(p)
  beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))
  y <- (X %*% beta_star + rnorm(n))
  
  g1 <- ApplyTriangleLassoHD( X = X, y = y, q = q, num_split = 1,
                              signal_index = signal_index, myseed = 1,mylambda=0.015)
  myresults=rbind(myresults,(c(s,g1$MDS_fdp,g1$MDS_power) ))
  print(myresults)
}
