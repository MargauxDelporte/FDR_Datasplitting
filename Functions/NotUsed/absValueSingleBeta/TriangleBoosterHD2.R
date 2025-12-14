#j=1
#model=lm
permR2TriangleBoostHD2<-function(data,j,model){
  dataPerm<-data[,-1]
  dataPerm[,j]<-sample(data[,j+1],replace=FALSE)
  names(dataPerm)=paste0('X',1:p)
  predictLM<-predict(model,newdata=data.matrix(dataPerm))
  rsquared=1-sum((data$y-predictLM)^2)/sum((data$y-mean(data$y))^2)
  return(rsquared)
}

ApplyTriangleBoostHD2<-function(X, y, q,myseed=1,myeta = 0.05,mymax_depth = 1,mylambda = 0.5,myalpha = 0.5, num_split=1,signal_index=signal_index){
  set.seed(myseed)
  amountTrain=0.333
  amountTest=1-amountTrain
  data<-data.frame(cbind(y,X))
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
  
  Xtrain=X[train_index,]
  names(Xtrain)=paste0('X',1:p)

  lm <- xgboost(
    data = data.matrix(Xtrain),
    label =y[train_index],
    booster = "gbtree",
    objective = "reg:squarederror",
    nrounds = 500,
    eta = myeta,
    max_depth = mymax_depth,
    lambda = mylambda,
    alpha = myalpha,
   # subsample = 0.7,
  #  colsample_bytree = 0.7,
    verbose = 0
  )
  summary(lm)
  remaining_index<-c(setdiff(c(1:n),train_index))
  sample_index1 <- sample(x = remaining_index, size = amountTest/2 * n, replace = F)
  sample_index2 <- setdiff(remaining_index, sample_index1)

  predict_TRAIN<-predict(lm,newdata=as.matrix(X[train_index,]))
  R2orig_TRAIN<-1-sum((y[train_index]-predict_TRAIN)^2)/sum((y[train_index]-mean(y[train_index]))^2)
  R2orig_TRAIN
  
  predictLM1<-predict(lm,newdata=as.matrix(X[sample_index1,]))
  predictLM2<-predict(lm,newdata=as.matrix(X[sample_index2,]))
  
  R2orig1<-1-sum((y[sample_index1]-predictLM1)^2)/sum((y[sample_index1]-mean(y[sample_index1]))^2)
  R2orig2<-1-sum((y[sample_index2]-predictLM2)^2)/sum((y[sample_index2]-mean(y[sample_index2]))^2)
  
  Rnew1<-sapply(1:ncol(X),function(j) permR2TriangleBoostHD2(data[sample_index1,],j,lm))
  Rnew2<-sapply(1:ncol(X),function(j) permR2TriangleBoostHD2(data[sample_index2,],j,lm))

  Diff1=R2orig1-Rnew1
  Diff2=R2orig2-Rnew2
  
  beta1=Diff1
  beta2=Diff2
  
  mirror<-sign(beta1*beta2)*(abs(beta1)+abs(beta2))
  hist(mirror[-signal_index])
  hist(mirror[signal_index])
  sort(signal_index)
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
  print(paste0('First R squared: ', round(R2orig1,3)))
  print(paste0('Second R squared: ', round(R2orig2,3)))
  print(paste0('DS_fdp = ', DS_fdp, ' DS_power = ', DS_power, ' MDS_fdp = ', MDS_fdp, ' MDS_power = ', MDS_power))
  return(list(DS_fdp = DS_fdp, DS_power = DS_power, MDS_fdp = MDS_fdp, MDS_power = MDS_power))
}

