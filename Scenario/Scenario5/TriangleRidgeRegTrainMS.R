library(MASS)
#j=1
permR2TriangleRidgeTrain<-function(mydata,j,model){
  dataPerm<-model.matrix(y ~ ., data = mydata)[, -1, drop = FALSE]
  dataPerm[,j]<-sample(mydata[,j+1],replace=FALSE)
  predictLM<-as.numeric(
    predict(
      model,
      newx = dataPerm,
      s = "lambda.min"
    )
  )
  rsquared=1-sum((mydata$y-predictLM)^2)/sum((mydata$y-mean(mydata$y))^2)
  return(rsquared)
}
myseed=5
ApplyTriangleRidgeTrain<-function(X, y, q,myseed,num_split=1,signal_index=signal_index){
  amountTrain=0.333
  amountTest=1-amountTrain
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
  
  # Build the predictor matrix and response vector
  xTrain <- model.matrix(y ~ ., data = dataTrain)[, -1, drop = FALSE]
  yTrain <- dataTrain$y
  
  # Fit ridge regression and select lambda by cross-validation
  ridge_fit <- cv.glmnet(
    x = xTrain,
    y = yTrain,
    family = "gaussian",
    type.measure = "mse",
  #  alpha = 0,
    nfolds = 25,
    standardize = F,
    relax=T
  )

  remaining_percent=1-amountTrain
  overlap=max(c(0,amountTest-remaining_percent))
  remaining_index<-c(setdiff(c(1:n),train_index),sample(train_index,size=overlap*n))
  sample_index1 <- sample(x = remaining_index, size = amountTest/2 * n, replace = F)
  sample_index2 <- setdiff(remaining_index, sample_index1)
  
  xTest1 <- model.matrix(y ~ ., data = data[sample_index1,])[, -1, drop = FALSE]
  xTest2 <- model.matrix(y ~ ., data = data[sample_index2,])[, -1, drop = FALSE]
  predictLM1 <- as.numeric(
    predict(
      ridge_fit,
      newx = xTest1,
      s = "lambda.min"
    )
  )
  predictLM2 <- as.numeric(
    predict(
      ridge_fit,
      newx = xTest2,
      s = "lambda.min"
    )
  )

  R2orig1<-1-sum((y[sample_index1]-predictLM1)^2)/sum((y[sample_index1]-mean(y[sample_index1]))^2)
  R2orig2<-1-sum((y[sample_index2]-predictLM2)^2)/sum((y[sample_index2]-mean(y[sample_index2]))^2)
  #print(c(R2orig1,R2orig2))
  Rnew1<-sapply(1:ncol(X),function(j) permR2TriangleRidgeTrain(mydata=data[sample_index1,],j,m=ridge_fit))
  Rnew2<-sapply(1:ncol(X),function(j) permR2TriangleRidgeTrain(mydata=data[sample_index2,],j,m=ridge_fit))
  
  beta1=R2orig1-Rnew1
  beta2=R2orig2-Rnew2
  

  mirror<-sign(beta1*beta2)*(abs(beta1)+abs(beta1))
  #hist(mirror)
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
  return(list(DS_fdp = DS_fdp, DS_power = DS_power, MDS_fdp = MDS_fdp, MDS_power = MDS_power))
}

