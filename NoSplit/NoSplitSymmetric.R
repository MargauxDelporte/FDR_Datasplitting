### High dimension linear model
rm(list = ls())

mywd='C:/Users/mde4023/OneDrive - Weill Cornell Medicine/0 Projects/FDR_Datasplitting'
mywd='C:/Users/mde4023/Documents/GitHub/FDR_Datasplitting'

setwd(mywd)
source('HelperFunctions.R')
source('TriangleLinRegTrainMS.R')
source('TriangleBoosterTrainMS.R')
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
#library(DSfdr)
library(ggplot2)
# library(igraph)
library(RColorBrewer)
library(openxlsx)
library(ggpubr)
library(viridis)

### algorithmic settings
num_split <- 1
n <-1500
p <- 250
p0 <- 25
q <- 0.1
#set.seed(124)(123) i=5
set.seed(456)
signal_index <- c()#sample(c(1:p), size = p0, replace = F)

permR2TriangleLinRegTrain<-function(data,j,model){
  dataPerm<-data[,-1]
  dataPerm[,j]<-sample(data[,j+1],replace=FALSE)
  predictLM<-predict(model,newdata=dataPerm)
  rsquared=1-sum((data$y-predictLM)^2)/sum((data$y-mean(data$y))^2)
  return(rsquared)
}
# my method that only returns the mirror statistic
ApplyLMNoSplit<-function(X, y, q,amountTrain=0.5,amountTest=1-amountTrain,myseed,num_split=1,signal_index=signal_index){
  set.seed(myseed)
  n <- dim(X)[1]; p <- dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  fdp <- rep(0, num_split)
  power <- rep(0, num_split)
  num_select <- rep(0, num_split)
  data<-data.frame(cbind(y,X))
  for(iter in 1:num_split){
    train_index<-sample(x = c(1:n), size = amountTrain * n, replace = F)
    dataTrain<-data[train_index,]
    colnames(dataTrain)<-c('y',paste0('X',1:p))
    colnames(data)<-c('y',paste0('X',1:p))
    
    lm<-lm(y~ ., data = dataTrain)
    remaining_percent=1-amountTrain
    overlap=max(c(0,amountTest-remaining_percent))
    remaining_index<-c(setdiff(c(1:n),train_index),sample(train_index,size=overlap*n))
    sample_index1 <- remaining_index
    
    
    predictLM1<-predict(lm,newdata=data.frame(data[sample_index1,]))
    
    R2orig1<-1-sum((y[sample_index1]-predictLM1)^2)/sum((y[sample_index1]-mean(y[sample_index1]))^2)
    
    Rnew1<-sapply(1:ncol(X),function(j) permR2TriangleLinRegTrain(data[sample_index1,],j,lm))
    
    Diff1=R2orig1-Rnew1
    
    sd_X1 <- apply(X[sample_index1, ], 2, sd)
    
    beta1=sign(Diff1)*sqrt(abs(Diff1))*sd(y)/sd_X1
    
    mirror<-beta1
    
    selected_index<-SelectFeatures(mirror,abs(mirror),0.1)
    
    ### number of selected variables j=1
    if(length(selected_index)!=0){
      num_select[iter] <- length(selected_index)
      inclusion_rate[iter, selected_index] <- 1/num_select[iter]
      
      ### calculate fdp and power
      result <- CalculateFDP_Power(selected_index, signal_index)
      fdp[iter] <- result$fdp
      power[iter] <- result$power
    }
  }
  
  ### single data-splitting (DS) result
  DS_fdp <- fdp[1]
  DS_power <- power[1]
  
  ### multiple data-splitting (MDS) result
  inclusion_rate <- apply(inclusion_rate, 2, mean)
  
  ### rank the features by the empirical inclusion rate
  feature_rank <- order(inclusion_rate)
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  if(length(feature_rank)!=0){
    null_feature <- numeric()
    
    ### backtracking 
    for(feature_index in 1:length(feature_rank)){
      if(sum(inclusion_rate[feature_rank[1:feature_index]]) > q){
        break
      }else{
        null_feature <- c(null_feature, feature_rank[feature_index])
      }
    }
    selected_index <- setdiff(feature_rank, null_feature)
    
    ### calculate fdp and power
    result <- CalculateFDP_Power(selected_index, signal_index)
    MDS_fdp <- result$fdp
    MDS_power <- result$power
  }
  else{
    MDS_fdp <- 0
    MDS_power <- 0
  }
  return(mirror)
  #list(DS_fdp = DS_fdp, DS_power = DS_power, MDS_fdp = MDS_fdp, MDS_power = MDS_power))
}

Compare_SignalStrenght=function(s){
  set.seed(s)
  ResultsDataFrame=data.frame()
  delta <- 5
  X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
  
  ### randomly generate the true beta i=4
  beta_star <- rep(0, p)
  beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
  
  ### generate y
  y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
  
  ###my own methods:
  g=ApplyLMNoSplit(X=as.data.frame(X), y, q=0.1,num_split=1, signal_index=signal_index, amountTrain=0.5, myseed = 1)
  
  ### save data
  return(g)}

Results=data.frame()
for(s in 1:25){
    Results=c(Results,Compare_SignalStrenght(s))
}

hist(as.numeric(Results))
