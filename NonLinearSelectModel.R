library(gbm)
library(ranger)
compute_R2 <- function(true, pred) {
  1 - sum((true - pred)^2) / sum((true - mean(true))^2)
}
Checkresults=function(lm){
  predict_TRAIN<-predict(lm,newdata=as.matrix(X[train_index,]))
  R2orig_TRAIN<-1-sum((y[train_index]-predict_TRAIN)^2)/sum((y[train_index]-mean(y[train_index]))^2)
  print(R2orig_TRAIN)
  
  predictLM1<-predict(lm,newdata=as.matrix(X[sample_index1,]))  
  R2orig1<-1-sum((y[sample_index1]-predictLM1)^2)/sum((y[sample_index1]-mean(y[sample_index1]))^2)
  print(R2orig1)
}


num_split <- 1
n <-1500
p <- 250
p0 <- 25
q <- 0.1
#set.seed(124)(123) i=5
set.seed(456)
signal_index <- sample(c(1:p), size = p0, replace = F)
ResultsDataFrame=data.frame()
delta <- 10
X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))

### randomly generate the true beta i=4
beta_star <- rep(0, p)
beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))*5

### generate y
y <- X^2%*%beta_star + rnorm(n, mean = 0, sd = 1)
y <-scale(y)


####create training and test set#####
amountTrain=0.333
amountTest=1-amountTrain
data<-data.frame(cbind(y,X))
names(data)=c('y',paste0('X',1:p))
train_index<-sample(x = c(1:n), size = amountTrain * n, replace = F)

dataTrain<-data[train_index,]
Xtrain=X[train_index,]
names(Xtrain)=paste0('X',1:p)

remaining_index<-c(setdiff(c(1:n),train_index))
sample_index1 <- sample(x = remaining_index, size = amountTest/2 * n, replace = F)
dataTest<-data[sample_index1,]


#fit the models


#xgboost

lm<-xgboost(data = data.matrix(Xtrain), label =y[train_index],nrounds=1000,lambda=5,eta=0.005,verbose=F,booster='gbtree')
Checkresults(lm)


#AdaBoost
gbm_mod <- gbm(y ~ .,
               data            = dataTrain,
               distribution    = "gaussian",
               n.trees         = 500,
               interaction.depth = 4,
               shrinkage       = 0.01,
               verbose         = FALSE)
pred_gbm_train <- predict(gbm_mod, newdata = dataTrain, n.trees = 500)
pred_gbm_test  <- predict(gbm_mod, newdata = dataTest,  n.trees = 500)
cat("GBM  R2 train:", compute_R2(dataTrain$y, pred_gbm_train), "\n")
cat("GBM  R2 test: ", compute_R2(dataTest$y,  pred_gbm_test),  "\n")


#ranger 

rf_sub <- ranger(
  formula         = y ~ ., 
  data            = dataTrain, 
  num.trees       = 500,
  mtry            = floor(sqrt(ncol(dataTrain) - 1)),  # typical default
  sample.fraction = 0.7,    # use 70% of the rows for each tree
  replace         = FALSE,  # without replacement for extra randomness
  importance      = "none"
)
pred_sub_train <- predict(rf_sub, data = dataTrain)$predictions
pred_sub_test <- predict(rf_sub, data = dataTest)$predictions
cat("GBM  R2 train:", compute_R2(dataTrain$y, pred_sub_train), "\n")
cat("GBM  R2 test: ", compute_R2(dataTest$y,  pred_sub_test),  "\n")
