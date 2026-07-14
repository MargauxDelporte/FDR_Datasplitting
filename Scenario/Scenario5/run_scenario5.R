###Linear model
rm(list = ls())

#mywd='/home/mde4023/FDR_Datasplitting'
mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
setwd(mywd)
source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Scenario/Scenario5/TriangleRidgeRegTrainMS.R'))
source(paste0(mywd,'/Functions Dai/knockoff.R'))
source(paste0(mywd,'/Functions Dai/analysis.R'))
source(paste0(mywd,'/Functions Dai/MBHq.R'))
source(paste0(mywd,'/Functions Dai/DS.R'))
source(paste0(mywd,'/Functions Dai/fdp_power.R'))
library(mvtnorm)
library(MASS)
library(glmnet)
library(knockoff)
library(mvtnorm)
library(parallel)
#install.packages("hdi") 
### algorithmic settings
num_split <- 10
n <- 800
p <- 100
p0 <- 10
q <- 0.05
rho <- 0.1
set.seed(456)
signal_index <- sample(c(1:p), size = p0, replace = F)

#i=8

Compare_SignalStrength=function(i,s){
  set.seed(s)
  delta <- i
  # simulate data
  n1 <- floor(n/2); n2 <- n - n1
  
  # Target equicorrelation matrix
  R0 <- matrix(rho, p, p)
  diag(R0) <- 1
  X1 <- rmvnorm(n1, rep(-1, p), sigma = R0)
  X2 <- rmvnorm(n2, rep(1, p), sigma = R0)
  X  <- rbind(X1, X2)
  X  <- rbind(X1, X2)
  beta_star <- rep(0, p)
  beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
  
  ### generate y
  y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
  
  ###my own methods:
  g=ApplyTriangleRidgeTrain(X=as.data.frame(X), y, q=q,num_split=num_split, signal_index=signal_index, myseed = 1)
  
  ResultsDataFrame=c('LinReg DS',i, as.numeric(g$DS_fdp),as.numeric(g$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('LinReg MS',i, as.numeric(g$MDS_fdp),as.numeric(g$MDS_power)))
  ResultsDataFrame
  return(ResultsDataFrame)
  ### Competition
  DS_result <- DS(X,y, num_split=num_split, q=q)
  ResultsDataFrame=rbind(ResultsDataFrame,c('DataSplitting',i,DS_result$DS_fdp,DS_result$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('MultipleDataSplitting',i,DS_result$MDS_fdp,DS_result$MDS_power))
  
  knockoff_result <- knockoff(X, y, q=q)
  ResultsDataFrame=rbind(ResultsDataFrame,c('Knockoff',i,knockoff_result$fdp,knockoff_result$power))
  
 # BH_result <- MBHq(X, y, q=q, num_split)
#  ResultsDataFrame=rbind(ResultsDataFrame,c('BH',i,BH_result$fdp,BH_result$power))
  
  ### save data
  return(ResultsDataFrame)}

Compare_SignalStrength(i=8,s=runif(1)*10000000)

# Load required packages
pkgs <- c('xgboost','gbm','ranger','MASS','glmnet','knockoff','mvtnorm','hdi',
          'foreach','doParallel')
lapply(pkgs, library, character.only = TRUE)

# === PARAMETER GRID ===
param_grid <- expand.grid(
  s = 1:200,
  i = seq(from = 7, to = 13, by = 1)
)
ncore=6
# === SET UP PARALLEL BACKEND ===
cl <- makeCluster(ncore-1)
# export working dir so workers can source
clusterExport(cl, 'mywd')
# have each worker source & load libraries
clusterEvalQ(cl, {
  setwd(mywd)
  source(file.path(mywd, '/Functions', 'HelperFunctions.R'))
  source(paste0(mywd,'/Scenario/Scenario5/TriangleLinRegTrainMS.R'))
  source(file.path(mywd, '/Functions Dai', 'knockoff.R'))
  source(file.path(mywd, '/Functions Dai', 'analysis.R'))
  source(file.path(mywd, '/Functions Dai', 'MBHq.R'))
  source(file.path(mywd, '/Functions Dai', 'DS.R'))
  source(file.path(mywd, '/Functions Dai', 'fdp_power.R'))
  lapply(c('xgboost','gbm','ranger','MASS','glmnet','knockoff','mvtnorm','hdi'),
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
  write.csv(chunk, file = paste0(mywd,"/Scenario/Scenario5/Temp2/",fname), row.names = FALSE)
  
  # return for final binding
  chunk
}
# === CLEANUP AND FINAL SAVE ===
stopCluster(cl)
warnings()
