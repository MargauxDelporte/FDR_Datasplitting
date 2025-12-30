### High dimension linear model
rm(list = ls())

mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
#mywd='C:/Users/mde4023/Documents/GitHub/FDR_Datasplitting'

setwd(mywd)

source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Functions/TriangleLassoHD.R'))
#source(paste0(mywd,'/Functions/TriangleElasticNetHD.R'))

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

### algorithmic settings
num_split <- 50
n <-400
p <- 500
p0 <- 25
q <- 0.1
set.seed(456)
signal_index <- sample(c(1:p), size = p0, replace = F)

#######set up the method for the comparison############# i=7 s=17 num_split=3
Compare_SignalStrength <- function(i, s) {
  set.seed(s)
  delta <- i
  # simulate data
  n1 <- floor(n/2); n2 <- n - n1
  X1 <- matrix(rnorm(n1*p, mean= -0.2), n1, p)
  X2 <- matrix(rnorm(n2*p, mean= 0.2), n2, p)
  X  <- rbind(X1, X2)
  beta_star <- numeric(p)
  beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))
  y <- (X %*% beta_star + rnorm(n))
  
  # run your custom methods
  g1 <- ApplyTriangleLassoHD( X = X, y = y, q = q, num_split = num_split,
                                 signal_index = signal_index, myseed = 1)
  print(g1)
  # FDR methods
  DS_result      <- DS(          X = X, y = y, q = q, num_split = num_split)
  knockoff_result<- knockoff(    X = X, y = y, q = q)
  BH_result      <- MBHq(        X = X, y = y, q = q, num_split = num_split)
  
  # init empty results df
  ResultsDataFrame <- data.frame(
    Method = character(),
    Delta  = numeric(),
    FDP    = numeric(),
    Power  = numeric(),
    stringsAsFactors = FALSE
  )
  
  # bind all rows
  ResultsDataFrame <- rbind(
    ResultsDataFrame,
    data.frame(Method = "Boost DS",                Delta = i, FDP = g1$DS_fdp,    Power = g1$DS_power),
    data.frame(Method = "Boost MS",                Delta = i, FDP = g1$MDS_fdp,   Power = g1$MDS_power),
    data.frame(Method = "DataSplitting",           Delta = i, FDP = DS_result$DS_fdp,  Power = DS_result$DS_power),
    data.frame(Method = "MultipleDataSplitting",   Delta = i, FDP = DS_result$MDS_fdp, Power = DS_result$MDS_power),
    data.frame(Method = "Knockoff",                Delta = i, FDP = knockoff_result$fdp, Power = knockoff_result$power),
    data.frame(Method = "Benjamini–Hochberg (BH)", Delta = i, FDP = BH_result$fdp,     Power = BH_result$power)
  )
  
  return(ResultsDataFrame)
}

#Compare_SignalStrength(8,8)

#######run the code#############
#Results=data.frame()
#for(s in 1:25){
#  for(i in seq(from=5,to=13,by=1)){
#  Results=rbind(Results,Compare_SignalStrength(i,s))
#  print(s)
#  }
#  print(Results)
#  }

library(parallel)

mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
setwd(mywd)

# Source helper and method files

source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Functions/TriangleLassoHD.R'))

# Dai’s routines
#source(file.path(mywd, 'Functions Dai', 'knockoff.R'))
source(file.path(mywd, 'Functions Dai', 'analysis.R'))
source(file.path(mywd, 'Functions Dai', 'MBHq.R'))
source(file.path(mywd, 'Functions Dai', 'DS.R'))
source(file.path(mywd, 'Functions Dai', 'fdp_power.R'))

# Load required packages
pkgs <- c('gbm','MASS','glmnet','knockoff','mvtnorm','hdi',
          'foreach','doParallel')
lapply(pkgs, library, character.only = TRUE)

# === PARAMETER GRID ===
param_grid <- expand.grid(
  s = 1:50,
  i = seq(from = 7, to = 13, by = 1)
)
ncore=detectCores()
# === SET UP PARALLEL BACKEND ===
cl <- makeCluster(ncore-1)
# export working dir so workers can source
clusterExport(cl, 'mywd')
# have each worker source & load libraries
clusterEvalQ(cl, {
  setwd(mywd)
  source(file.path(mywd, 'Functions', 'TriangleLassoHD.R'))
  source(file.path(mywd, 'Functions', 'HelperFunctions.R'))
  source(file.path(mywd, 'Functions Dai', 'knockoff.R'))
  source(file.path(mywd, 'Functions Dai', 'analysis.R'))
  source(file.path(mywd, 'Functions Dai', 'MBHq.R'))
  source(file.path(mywd, 'Functions Dai', 'DS.R'))
  source(file.path(mywd, 'Functions Dai', 'fdp_power.R'))
  lapply(c('gbm','MASS','glmnet','knockoff','mvtnorm','hdi'),
         library, character.only = TRUE)
})
registerDoParallel(cl)

# === RUN IN PARALLEL AND WRITE OUT ===
results_list <- foreach(
  k = seq_len(nrow(param_grid)),
  .packages = pkgs,
  .combine  = rbind
) %dopar% {
  s_val <- param_grid$s[k]
  i_val <- param_grid$i[k]
  
  # compute chunk of results
  chunk <- Compare_SignalStrength(i = i_val, s = s_val)
  
  # write out this chunk immediately
  fname <- sprintf("Results_s%02d_i%02d.csv", s_val, i_val)
  write.csv(chunk, file = paste0("C:/Users/mde4023/Downloads/FDR_Datasplitting/Results/ResultsHD scenario/Temp2/",fname), row.names = FALSE)
  
  # return for final binding
  chunk
}

# === CLEANUP AND FINAL SAVE ===
stopCluster(cl)
warnings()
# combine all and save full dataset
Results <- results_list
write.csv(Results, file = "ResultsHDScenario.csv", row.names = FALSE)

