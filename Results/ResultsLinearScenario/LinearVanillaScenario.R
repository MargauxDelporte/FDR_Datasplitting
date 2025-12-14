### High dimension linear model
rm(list = ls())

mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
#mywd='C:/Users/mde4023/Documents/GitHub/FDR_Datasplitting'
#C:\Users\mde4023\Downloads\FDR_Datasplitting

setwd(mywd)
source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Functions/TriangleLinRegTrainMS.R'))
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
library(dplyr)
library(ggplot2)
library(ggpubr)


### algorithmic settings
num_split <- 50
n <-1500
p <- 250
p0 <- 25
q <- 0.1
#set.seed(124)(123) i=5
set.seed(456)
signal_index <- sample(c(1:p), size = p0, replace = F)

################################Triange: train test test each 30####################################
Compare_SignalStrength=function(i,s){
  set.seed(s)
  ResultsDataFrame=data.frame()
  delta <- i
  X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
  
  ### randomly generate the true beta i=4
  beta_star <- rep(0, p)
  beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
  
  ### generate y
  y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
  
  ###my own methods:
  g=ApplyTriangleLinRegTrain(X=as.data.frame(X), y, q=0.1,num_split=num_split, signal_index=signal_index, amountTrain=0.333, myseed = 1)
  
  ResultsDataFrame=c('LinReg DS',i, as.numeric(g$DS_fdp),as.numeric(g$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('LinReg MS',i, as.numeric(g$MDS_fdp),as.numeric(g$MDS_power)))
  
  ### Competition
  DS_result <- DS(X,y, num_split=num_split, q=0.1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('DataSplitting',i,DS_result$DS_fdp,DS_result$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('MultipleDataSplitting',i,DS_result$MDS_fdp,DS_result$MDS_power))
  
  knockoff_result <- knockoff(X, y, q=0.1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('Knockoff',i,knockoff_result$fdp,knockoff_result$power))
  
  BH_result <- MBHq(X, y, q=0.1, num_split)
  ResultsDataFrame=rbind(ResultsDataFrame,c('BH',i,BH_result$fdp,BH_result$power))
  
  ### save data
  return(ResultsDataFrame)}

library(parallel)
mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
setwd(mywd)

# Source helper and method files
source(file.path(mywd, 'Functions','HelperFunctions.R'))
source(file.path(mywd, 'Functions', 'TriangleLinRegTrainMS.R'))

# Dai’s routines
source(file.path(mywd, 'Functions Dai', 'knockoff.R'))
source(file.path(mywd, 'Functions Dai', 'analysis.R'))
source(file.path(mywd, 'Functions Dai', 'MBHq.R'))
source(file.path(mywd, 'Functions Dai', 'DS.R'))
source(file.path(mywd, 'Functions Dai', 'fdp_power.R'))

# Load required packages
pkgs <- c('xgboost','gbm','ranger','MASS','glmnet','knockoff','mvtnorm','hdi',
          'foreach','doParallel')
lapply(pkgs, library, character.only = TRUE)

# === PARAMETER GRID ===
param_grid <- expand.grid(
  s = 1:50,
  i = seq(from = 9, to = 13, by = 1)
)
ncore=detectCores()
# === SET UP PARALLEL BACKEND ===
cl <- makeCluster(ncore-1)
# export working dir so workers can source
clusterExport(cl, 'mywd')
# have each worker source & load libraries
clusterEvalQ(cl, {
  setwd(mywd)
  source(file.path(mywd, 'Functions', 'HelperFunctions.R'))
  source(file.path(mywd, 'Functions', 'TriangleLinRegTrainMS.R'))
  source(file.path(mywd, 'Functions Dai', 'knockoff.R'))
  source(file.path(mywd, 'Functions Dai', 'analysis.R'))
  source(file.path(mywd, 'Functions Dai', 'MBHq.R'))
  source(file.path(mywd, 'Functions Dai', 'DS.R'))
  source(file.path(mywd, 'Functions Dai', 'fdp_power.R'))
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
  write.csv(chunk, file = paste0(mywd,"/Results/ResultsLinearScenario/Temp2/",fname), row.names = FALSE)
  
  # return for final binding
  chunk
}
# === CLEANUP AND FINAL SAVE ===
stopCluster(cl)
warnings()
