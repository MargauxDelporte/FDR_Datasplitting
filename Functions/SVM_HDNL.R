#j=1
#model=lm
permR2SVM<-function(data,j,model){
  dataPerm<-data[,-1]
  dataPerm[,j]<-sample(data[,j+1],replace=FALSE)
  names(dataPerm)=paste0('X',1:p)
  predictLM<-predict(model,newdata=(dataPerm))
  rsquared=1-sum((data$y-predictLM)^2)/sum((data$y-mean(data$y))^2)
  return(rsquared)
}

ApplySVM<-function(X, y, q,myseed=1,myparms,num_split=1,signal_index=signal_index){
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
    
    Xtrain=dataTrain[,-1]
    ytrain=y[train_index]
    fit=krr_model <- ksvm(
      x = as.matrix(Xtrain),
      y = ytrain,
      kernel = polydot(degree = myparms$degree, scale = myparms$scale, offset = myparms$offset), lambda = myparms$lambda
    )
    
    remaining_index<-c(setdiff(c(1:n),train_index))
    sample_index1 <- sample(x = remaining_index, size = amountTest/2 * n, replace = F)
    sample_index2 <- setdiff(remaining_index, sample_index1)
    
    predict_TRAIN<-predict(fit,newdata=Xtrain)
    R2orig_TRAIN<-1-sum((y[train_index]-predict_TRAIN)^2)/sum((y[train_index]-mean(y[train_index]))^2)
    R2orig_TRAIN
    
    dataTest1<-data[sample_index1,-1]
    dataTest2<-data[sample_index2,-1]
    predictLM1<-predict(fit,newdata=dataTest1)
    predictLM2<-predict(fit,newdata=dataTest2)
    
    R2orig1<-1-sum((y[sample_index1]-predictLM1)^2)/sum((y[sample_index1]-mean(y[sample_index1]))^2)
    R2orig2<-1-sum((y[sample_index2]-predictLM2)^2)/sum((y[sample_index2]-mean(y[sample_index2]))^2)
    
    Rnew1<-sapply(1:ncol(X),function(j) permR2SVM(data[sample_index1,],j,fit))
    Rnew2<-sapply(1:ncol(X),function(j) permR2SVM(data[sample_index2,],j,fit))
    
    Diff1=R2orig1-Rnew1
    Diff2=R2orig2-Rnew2
    
    beta1=Diff1
    beta2=Diff2
    
    mirror<-sign(beta1*beta2)*(abs(beta1)+abs(beta2))
    hist(mirror[-signal_index])
    hist(mirror[signal_index])
    selected_index<-SelectFeatures(mirror,abs(mirror),q)
    
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
  print(hist(mirror[-signal_index]))
  print(paste0('First R squared: ', round(R2orig1,3)))
  print(paste0('Second R squared: ', round(R2orig2,3)))
  print(paste0('DS_fdp = ', DS_fdp, ' DS_power = ', DS_power, ' MDS_fdp = ', MDS_fdp, ' MDS_power = ', MDS_power))
  return(list(DS_fdp = DS_fdp, DS_power = DS_power, MDS_fdp = MDS_fdp, MDS_power = MDS_power))
}


