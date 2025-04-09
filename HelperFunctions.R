### Supporting functions
SelectFeatures <- function(mm, ww, q){
  ### mm: mirror statistics mm=M
  ### ww: absolute value of mirror statistics ww=abs(M)
  ### q:  FDR control level q=0.1
  cutoff_set <- max(ww)
  for(t in ww){
    ps <- length(mm[mm > t])
    ng <- length(na.omit(mm[mm < -t]))
    rto <- (ng + 1)/max(ps, 1)
    if(rto <= q){
      cutoff_set <- c(cutoff_set, t)
    }
  }
  cutoff <- min(cutoff_set)
  selected_index <- which(mm > cutoff)
  
  return(selected_index)
}

CalculateFDP_Power <- function(selected_index, signal_index){
  num_selected <- length(selected_index)
  tp <- length(intersect(selected_index, signal_index))
  fp <- num_selected - tp
  fdp <- fp / max(num_selected, 1)
  power <- tp / length(signal_index)
  return(list(fdp = fdp, power = power))
}

ChangePredictFeature=function(X,ncovar){
  Predictor <- as.matrix(X[,-ncovar])
  y=as.vector(X[,ncovar])
  Tempdata=as.data.frame(cbind(y,Predictor))
  X[,ncovar]=rfsrc(y~.,data=Tempdata)$predicted
  X=as.data.frame(X)
  names(X)=paste0('X',1:p)
  return(X)
}

###ncovar=10
ChangePredictFeature=function(X,ncovar){
  Predictor <- as.matrix(X[,-ncovar])
  y=as.vector(X[,ncovar])
  Tempdata=as.data.frame(cbind(y,Predictor))
  X[,ncovar]=rfsrc(y~.,data=Tempdata)$predicted
  X=as.data.frame(X)
  names(X)=paste0('X',1:p)
  return(X)
}

CreateDatasetPredictedFeatureBoost=function(X){
  p=ncol(X)
  AltDataset=data.frame(rep(0,nrow(X)))
  for(ncovar in 1:p){
    y=as.vector(X[,ncovar])
    TempModel=xgboost(data = data.matrix(X[,-ncovar]), label =y,nrounds=5)
    AltDataset=cbind(AltDataset,predict(TempModel,as.matrix(X[,-ncovar])))
    print(ncovar/p)
  }
  AltDataset=AltDataset[,-1]
  names(AltDataset)=paste0('X',1:p)
  return(AltDataset)
}
CreateDatasetPermutedfeatures=function(X){
  p=ncol(X)
  AltDataset=data.frame(rep(0,nrow(X)))
  for(ncovar in 1:p){
    AltDataset=cbind(AltDataset,sample(X[,ncovar]))
  }
  AltDataset=AltDataset[,-1]
  names(AltDataset)=paste0('X',1:p)
  return(AltDataset)
}
CreateDatasetPredictedFeaturesLasso=function(X){
  p=ncol(X)
  AltDataset=data.frame(rep(0,nrow(X)))
  for(ncovar in 1:p){
    y=as.vector(X[,ncovar])
    x=data.matrix(X[,-ncovar])
    cvfit <- cv.glmnet(data.matrix(x), y, type.measure = "mse", nfolds = 10)
    lambda <- cvfit$lambda.min
    TempModel=glmnet(data.matrix(x), y, family = "gaussian", alpha = 1, lambda = lambda)
    AltDataset=cbind(AltDataset,predict(TempModel,as.matrix(X[,-ncovar])))
    print(ncovar/p)
  }
  AltDataset=AltDataset[,-1]
  names(AltDataset)=paste0('X',1:p)
  return(AltDataset)
}
CreateDatasetPredictedFeaturesForest=function(X){
  p=ncol(X)
  AltDataset=data.frame(rep(0,nrow(X)))
  for(ncovar in 1:p){
    y=as.vector(X[,ncovar])
    x=data.matrix(X[,-ncovar])
    TempData=as.data.frame(cbind(y,x))
    TempModel=rfsrc(y ~ .,data=TempData)
    AltDataset=cbind(AltDataset,predict(TempModel,TempData[,-1])$predicted)
    print(ncovar/p)
  }
  AltDataset=AltDataset[,-1]
  names(AltDataset)=paste0('X',1:p)
  return(AltDataset)
}

