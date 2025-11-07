# Load libraries
library(nnet)
library(caret)
library(MASS)
library(glmnet)
library(knockoff)
library(mvtnorm)
library(hdi)
library(randomForest)
library(dplyr)
library(psych) #Partial correlation
library(xgboost)
library(ggplot2)
library(igraph)
library(RColorBrewer)
library(openxlsx)
library(ggpubr)
library(viridis)
test=function(){
  predictions1 <- mymodel %>% predict(train_x)
  SSE1 <- sum((train_y - predictions1)^2)
  TSE1 <- sum((train_y - mean(train_y))^2)
  R2orig1 <- 1 - SSE1 / TSE1
  R2orig1
  
  predictions1 <- mymodel %>% predict(test_x)
  SSE1 <- sum((test_y - predictions1)^2)
  TSE1 <- sum((test_y - mean(test_y))^2)
  R2orig12 <- 1 - SSE1 / TSE1
  R2orig12
  return(c(R2orig1,R2orig12))
}

setwd('C:/Users/mde4023/Documents/GitHub/Mirror_DS_MSE')
source('HelperFunctions.R')
source('C:/Users/mde4023/Documents/GitHub/Mirror_DS_MSE/TriangleBoosterTrainMS.R')
source('C:/Users/mde4023/Documents/GitHub/Mirror_DS_MSE/Functions Dai/DS.R')
source('C:/Users/mde4023/Documents/GitHub/Mirror_DS_MSE/Functions Dai/analysis.R')
source('C:/Users/mde4023/Documents/GitHub/Mirror_DS_MSE/Functions Dai/MBHq.R')
source('C:/Users/mde4023/Documents/GitHub/Mirror_DS_MSE/Functions Dai/knockoff.R')
source('C:/Users/mde4023/Documents/GitHub/Mirror_DS_MSE/Functions Dai/fdp_power.R')



###################################################
# 1. Data Generation
###################################################
set.seed(123)
n     <- 5000
p     <- 250
p0    <- 20
delta=10

# Generate X ~ N(0, I)
X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))

# Suppose signal_index are the first p0 for demonstration
signal_index <- 1:p0

# True quadratic coefficients
alpha_star <- rep(0, p)
alpha_star[signal_index] <- rnorm(p0, mean = 0, sd = delta * sqrt(log(p) / n))

# Create elementwise squares of X
X2 <- X^2

# Generate y with both linear and quadratic terms, then scale
y <- X2 %*% alpha_star + rnorm(n, mean = 0, sd = 1)#+X %*% beta_star + 
y <- scale(y)  # Standardize response

###################################################
# 2. Combine X and X^2 into a single data frame
###################################################
# Rename columns to avoid duplicates
colnames(X)  <- paste0("X", 1:p)

# Optionally, we can scale each feature to zero mean, unit variance
X  <- scale(X)

# Combine into one data frame
df <- data.frame(y = as.numeric(y), X)
###################################################
# 3. Split into Training and Test sets
###################################################

#repeat this 150 times
Results=data.frame()
SimulationFunction=function(s){
  set.seed(s)
  ResultsDataFrame=data.frame()
  X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
  delta=10
  # Suppose signal_index are the first p0 for demonstration
  signal_index <- 1:p0
  
  # True linear coefficients
  beta_star <- rep(0, p)
  beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta * sqrt(log(p) / n))
  
  # True quadratic coefficients
  alpha_star <- rep(0, p)
  alpha_star[signal_index] <- rnorm(p0, mean = 0, sd = delta * sqrt(log(p) / n))
  
  # Create elementwise squares of X
  X2 <- X^2
  
  # Generate y with both linear and quadratic terms, then scale
  y <- X %*% beta_star + X2 %*% alpha_star + rnorm(n, mean = 0, sd = 1)
  y <- scale(y)  # Standardize response
  
  # Rename columns to avoid duplicates
  colnames(X)  <- paste0("X", 1:p)
  
  # Apply function
  g2=ApplyTriangleBoostTrain(X, y, q=0.1,num_split=10,mybooster='gbtree', signal_index=signal_index, amountTrain=0.333, myseed = 1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('GBtree DS',g2$DS_fdp,g2$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('GBtree MS',g2$MDS_fdp,g2$MDS_power))
  
  ### Competition
  DS_result <- DS(X,y, num_split=10, q=0.1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('DataSplitting',DS_result$DS_fdp,DS_result$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('MultipleDataSplitting',DS_result$MDS_fdp,DS_result$MDS_power))
  
  knockoff_result <- knockoff(X, y, q=0.1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('Knockoff',knockoff_result$fdp,knockoff_result$power))
  
  BH_result <- MBHq(X, y, q=0.1,num_split=10)
  ResultsDataFrame=rbind(ResultsDataFrame,c('BH',BH_result$fdp,BH_result$power))
  
  names(ResultsDataFrame)=c('Method', 'FDR','Power')
  ### save data
  return(ResultsDataFrame)}

Results=data.frame()
for(s in 1:150){
  Results=rbind(Results,SimulationFunction(s))
  print(Results)
  print(s)
}
getwd()
library(writexls)
write.csv(Results,file='ResultsNonLinearTriangle.csv')

Results2=Results
head(Results2)
Results2$FDR =round(as.numeric(Results2$FDR),3)
Results2$Power=round(as.numeric(Results2$Power),2)

resultsagg <- Results2 %>%
  group_by(Method) %>%
  summarize(
    Avg_FDR = mean(FDR),
    Avg_Power = mean(Power),
  )
resultsagg


########################################
#######FOR KERAS########################
########################################
