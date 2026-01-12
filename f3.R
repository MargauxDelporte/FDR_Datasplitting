### High dimension linear model
rm(list = ls())

mywd='/home/mde4023/FDR_Datasplitting'
setwd(mywd)

source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Scenario/Scenario3b/TriangleLassoHD.R'))

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
num_split <- 50
n <-400
p <- 500
p0 <- 25
q <- 0.1
set.seed(456)
signal_index <- sample(c(1:p), size = p0, replace = F)

  set.seed(5)
  delta <- 10
  # simulate data
  n1 <- floor(n/2); n2 <- n - n1
  X1 <- matrix(rnorm(n1*p, mean= -1), n1, p)
  X2 <- matrix(rnorm(n2*p, mean= 1), n2, p)
  X  <- rbind(X1, X2)
  beta_star <- numeric(p)
  beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))
  y <- (X %*% beta_star + rnorm(n))
  
  set.seed(123)
  amountTrain=0.333
  amountTest=1-amountTrain
  data<-data.frame(cbind(y,X))
  n <- dim(X)[1]; p <- dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  fdp <- rep(0, num_split)
  power <- rep(0, num_split)
  num_select <- rep(0, num_split)
  data<-data.frame(cbind(y,X))

    train_index<-sample(x = c(1:n), size = amountTrain * n, replace = F)
    dataTrain<-data[train_index,]
    colnames(dataTrain)<-c('y',paste0('X',1:p))
    colnames(data)<-c('y',paste0('X',1:p))
    
    Xtrain=X[train_index,]
    names(Xtrain)=paste0('X',1:p)
    
    # crossvalidation lasso
    cv_lasso <- cv.glmnet(x = as.matrix(Xtrain), y = y[train_index], alpha = 1)
    
    best_lambda <- cv_lasso$lambda.1se*1.5
    # Fit final model with best lambda
    lm <- glmnet(x = as.matrix(Xtrain), y = y[train_index], alpha = 1, lambda = best_lambda)
    
    remaining_index<-c(setdiff(c(1:n),train_index))
    sample_index1 <- sample(x = remaining_index, size = amountTest/2 * n, replace = F)
    sample_index2 <- setdiff(remaining_index, sample_index1)
    
    predict_TRAIN<-predict(lm,newx=as.matrix(X[train_index,]))
    R2orig_TRAIN<-1-sum((y[train_index]-predict_TRAIN)^2)/sum((y[train_index]-mean(y[train_index]))^2)
    R2orig_TRAIN
    
    predictLM1<-predict(lm,newx=as.matrix(X[sample_index1,]))
    predictLM2<-predict(lm,newx=as.matrix(X[sample_index2,]))
    
    R2orig1<-1-sum((y[sample_index1]-predictLM1)^2)/sum((y[sample_index1]-mean(y[sample_index1]))^2)
    R2orig2<-1-sum((y[sample_index2]-predictLM2)^2)/sum((y[sample_index2]-mean(y[sample_index2]))^2)
    
    Rnew1<-sapply(1:ncol(X),function(j) permR2TriangleLassHD(data[sample_index1,],j,lm))
    Rnew2<-sapply(1:ncol(X),function(j) permR2TriangleLassHD(data[sample_index2,],j,lm))
    
    beta1=R2orig1-Rnew1
    beta2=R2orig2-Rnew2
    
    
    mirror<-sign(beta1*beta2)*(abs(beta1)+abs(beta1))
    hist(mirror)
    selected_index<-SelectFeatures(mirror,abs(mirror),0.1)
    
    ### number of selected variables j=1
    R2orig1; R2orig2
      result <- CalculateFDP_Power(selected_index, signal_index)
      fdp <- result$fdp
      power <- result$power
      fdp
      power

