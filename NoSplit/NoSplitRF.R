library(MASS)
library(randomForest)

#mydata=data[remaining_index,]  j=1 model=rf
permR2RFNoSplit<-function(mydata,j,model){
  dataPerm<-mydata[,-1]
  dataPerm[,j]<-sample(mydata[,j+1],replace=FALSE)
  names(dataPerm)=c()
  predictLM<-predict(model, newdata = as.matrix(dataPerm))
  rsquared=1-sum((mydata$y-predictLM)^2)/sum((mydata$y-mean(mydata$y))^2)
  return(rsquared)
}

ApplyRFNoSplit<-function(X, y, q,amountTrain=0.5,myseed,num_split=1,signal_index=signal_index){
  set.seed(myseed)
  n <- dim(X)[1]; p <- dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  fdp <- rep(0, num_split)
  power <- rep(0, num_split)
  num_select <- rep(0, num_split)
  data<-data.frame(cbind(y,X))
  names(data)=c('y',paste0('X',1:p))
  for(iter in 1:num_split){
  train_index<-sample(x = c(1:n), size = amountTrain * n, replace = F)
  dataTrain<-data[train_index,]
  colnames(dataTrain)<-c('y',paste0('X',1:p))
  colnames(data)<-c('y',paste0('X',1:p))
  names(X)=c(paste0('X',1:p))
  rf<-randomForest(x = as.matrix(X[train_index, ]), y = y[train_index])
  remaining_percent=1-amountTrain
  remaining_index<-c(setdiff(c(1:n),train_index),sample(train_index,size=remaining_percent*n))

  predictLM1<-predict(rf, newdata = as.matrix(X[remaining_index, ]))

  R2orig1<-1-sum((y[remaining_index]-predictLM1)^2)/sum((y[remaining_index]-mean(y[remaining_index]))^2)

  Rnew1<-sapply(1:ncol(X),function(j) permR2RFNoSplit(data[remaining_index,],j,rf))

  Diff1=R2orig1-Rnew1

  sd_X1 <- apply(X[sample_index1, ], 2, sd)

  beta1=sign(Diff1)*sqrt(abs(Diff1))*sd(y)/sd_X1

  mirror<-beta1
  #mirror2=c()
  #mirror2=c(mirror2,mirror)
  #hist(mirror2)
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
  return(list(DS_fdp = DS_fdp, DS_power = DS_power, MDS_fdp = MDS_fdp, MDS_power = MDS_power))
}

