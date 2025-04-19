library(reticulate)
use_condaenv("r-tf-py310", required = TRUE)
library(keras)
library(tensorflow)

model <- keras_model_sequential() %>%
  # 1) square every input feature
  layer_lambda(
    f         = function(x) k_square(x),
    input_shape = p,
    name      = "quadratic_layer"
  ) %>%
  # 2) linear regression on the squared features
  layer_dense(
    units      = 1,
    activation = "linear",
    name       = "output"
  )
permR2TriangleNNTrain<-function(data,j,model){
  dataPerm<-data[,-1]
  dataPerm[,j]<-sample(data[,j+1],replace=FALSE)
  names(dataPerm)=paste0('X',1:p)
  predictLM<-compute(nn, dataPerm)$net.result
  rsquared=1-sum((data$y-predictLM)^2)/sum((data$y-mean(data$y))^2)
  return(rsquared)
}

ApplyTrianglNNTrain<-function(X, y, q,myseed,num_split=1,signal_index=signal_index){
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
  feature_names <- paste0('X',1:p)
  fmla <- as.formula(paste("y ~", paste(feature_names, collapse = " + ")))
  
  # 5) Fit the neural network
  scaled_train <- as.data.frame(scale(dataTrain))
  scaled_train$y <- scale(dataTrain$y)
  nn_fixed <- neuralnet(
    formula       = fmla,
    data          = scaled_train,
    hidden        = 50,              # much smaller than 250
    act.fct       = function(x) x^2, # quadratic
    linear.output = TRUE,
    
    algorithm     = "backprop",      # faster, more stable
    learningrate  = 0.005,           # moderate rate
    threshold     = 1e-4,            # stop when SSE change < 1e-4
    stepmax       = 1e6,             # allow up to a million steps
    rep           = 3,               # integer: run 3 random starts
    lifesign      = "minimal"
  )
  nn = neuralnet(
    fmla,
    data=dataTrain,
    hidden=c(240),
    act.fct       = function(x) x^2,
    linear.output = FALSE
  )
  
  remaining_index<-c(setdiff(c(1:n),train_index))
  sample_index1 <- sample(x = remaining_index, size = amountTest/2 * n, replace = F)
  sample_index2 <- setdiff(remaining_index, sample_index1)

  predictLM1 <- compute(nn, X[sample_index1,])$net.result
  predictLM2 <- compute(nn, X[sample_index2,])$net.result
  
  R2orig1<-1-sum((y[sample_index1]-predictLM1)^2)/sum((y[sample_index1]-mean(y[sample_index1]))^2)
  R2orig2<-1-sum((y[sample_index2]-predictLM2)^2)/sum((y[sample_index2]-mean(y[sample_index2]))^2)
  
  Rnew1<-sapply(1:ncol(X),function(j) permR2TriangleNNTrain(data[sample_index1,],j,lm))
  Rnew2<-sapply(1:ncol(X),function(j) permR2TriangleNNTrain(data[sample_index2,],j,lm))

  Diff1=R2orig1-Rnew1
  Diff2=R2orig2-Rnew2
  
  sd_X1 <- apply(X[sample_index1, ], 2, sd)
  sd_X2 <- apply(X[sample_index2, ], 2, sd)
  
  
  beta1=sign(Diff1)*sqrt(abs(Diff1))*sd(y)/sd_X1
  beta2=sign(Diff2)*sqrt(abs(Diff2))*sd(y)/sd_X2
  
  mirror<-sign(beta1*beta2)*(abs(beta1))
  hist(mirror[-signal_index])
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
  print(paste0('First R squared: ', round(Rnew1[1],3)))
  print(paste0('Second R squared: ', round(Rnew2[1],3)))
  print(paste0('DS_fdp = ', DS_fdp, ' DS_power = ', DS_power, ' MDS_fdp = ', MDS_fdp, ' MDS_power = ', MDS_power))
  return(list(DS_fdp = DS_fdp, DS_power = DS_power, MDS_fdp = MDS_fdp, MDS_power = MDS_power))
}

