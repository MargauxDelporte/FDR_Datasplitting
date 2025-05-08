library(caret)
library(xgboost)
library(lightgbm)
library(MASS)
### algorithmic settings
num_split <- 1
n <-1500
p <- 2000
p0 <- 25
q <- 0.1

#sample the data
set.seed(456)
signal_index <- sample(c(1:p), size = p0, replace = F)
delta <- 10
X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
beta_star <- numeric(p)
beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))*100
y <- scale(X^2 %*% beta_star + rnorm(n))



#regular start method
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
remaining_index<-c(setdiff(c(1:n),train_index))
sample_index1 <- sample(x = remaining_index, size = amountTest/2 * n, replace = F)

dtrain <- xgb.DMatrix(data = as.matrix(Xtrain), label = y[train_index])
dtest <- xgb.DMatrix(data = as.matrix(X[sample_index1, ]))
y_true <- y[sample_index1]

# define grid
#param_grid <- expand.grid(
# myeta=0.2,mymax_depth = 1,mylambda = 0.1,myalpha=0.1
train_data <- lgb.Dataset(data = data.matrix(Xtrain), label = y[train_index])

# Set parameters
param_grid <- expand.grid(
  learning_rate = c(0.01, 0.05, 0.1),
  max_depth = c(3, 5, 7),
  num_leaves = c(15, 31, 63),  # must be < 2^max_depth
  lambda_l1 = c(0, 1, 5),
  lambda_l2 = c(0, 1, 5),
  min_data_in_leaf = c(10, 20, 50),
  feature_fraction = c(0.6, 0.8, 1.0),
  bagging_fraction = c(0.6, 0.8, 1.0),
  bagging_freq = c(1)  # how often to apply bagging
)
as.list(param_grid[2, ])
results <- data.frame()


# Run grid search i=1
for (i in 1:nrow(param_grid)) {
  # Convert row to named list and add required fields
  params <- as.list(param_grid[i, ])
  params$objective <- "regression"
  params$metric <- "l2"
  params$verbosity <- -1  # silent mode
  
  model <- lgb.train(
    params = params,
    data = train_data,
    nrounds = 500,
    verbose = -1
  )
  
  preds <- predict(model, data.matrix(X[sample_index1, ]))
  R2 <- 1 - sum((y_true - preds)^2) / sum((y_true - mean(y_true))^2)
  newrow=c(round(i), R2)
  results <- rbind(results, newrow)
  print(newrow)
}

results <- results[order(-results[,2]), ]
print(head(results))  


permR2TriangleBoostHD2<-function(data,j,model){
  dataPerm<-data[,-1]
  dataPerm[,j]<-sample(data[,j+1],replace=FALSE)
  names(dataPerm)=paste0('X',1:p)
  predictLM<-predict(model, data.matrix(dataPerm))
  rsquared=1-sum((data$y-predictLM)^2)/sum((data$y-mean(data$y))^2)
  return(rsquared)
}

params <- as.list(param_grid[2, ])
params$objective <- "regression"
params$metric <- "l2"
params$verbosity <- -1  # silent mode

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
    train_data <- lgb.Dataset(data = data.matrix(Xtrain), label = y[train_index])
    
    lm <- lgb.train(
      params = params,
      data = train_data,
      nrounds = 500,
      verbose = -1
    )
    
    remaining_index<-c(setdiff(c(1:n),train_index))
    sample_index1 <- sample(x = remaining_index, size = amountTest/2 * n, replace = F)
    sample_index2 <- setdiff(remaining_index, sample_index1)
    
    predict_TRAIN<-predict(lm,as.matrix(X[train_index,]))
    R2orig_TRAIN<-1-sum((y[train_index]-predict_TRAIN)^2)/sum((y[train_index]-mean(y[train_index]))^2)
    R2orig_TRAIN
    
    predictLM1<-predict(lm,data.matrix(X[sample_index1,]))
    predictLM2<-predict(lm,data.matrix(X[sample_index2,]))
    
    R2orig1<-1-sum((y[sample_index1]-predictLM1)^2)/sum((y[sample_index1]-mean(y[sample_index1]))^2)
    R2orig2<-1-sum((y[sample_index2]-predictLM2)^2)/sum((y[sample_index2]-mean(y[sample_index2]))^2)
    
    Rnew1<-sapply(1:ncol(X),function(j) permR2TriangleBoostHD2(data[sample_index1,],j,lm))
    Rnew2<-sapply(1:ncol(X),function(j) permR2TriangleBoostHD2(data[sample_index2,],j,lm))
    
    Diff1=R2orig1-Rnew1
    Diff2=R2orig2-Rnew2
    
    sd_X1 <- apply(X[sample_index1, ], 2, sd)
    sd_X2 <- apply(X[sample_index2, ], 2, sd)
    
    
    beta1=sign(Diff1)*sqrt(abs(Diff1))*sd(y)/sd_X1
    beta2=sign(Diff2)*sqrt(abs(Diff2))*sd(y)/sd_X2
    
    mirror<-sign(beta1*beta2)*(abs(beta1))
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
  #print(paste0('First R squared: ', round(R2orig1,3)))
  #print(paste0('Second R squared: ', round(R2orig2,3)))
  # print(paste0('DS_fdp = ', DS_fdp, ' DS_power = ', DS_power, ' MDS_fdp = ', MDS_fdp, ' MDS_power = ', MDS_power))
  return(list(DS_fdp = DS_fdp, DS_power = DS_power, MDS_fdp = MDS_fdp, MDS_power = MDS_power))
}


ApplyTriangleBoostHD2(X=X, y=y, q=q, myseed=1, num_split=1, signal_index=signal_index)
