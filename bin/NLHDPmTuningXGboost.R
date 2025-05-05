library(caret)
library(xgboost)

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
param_grid <- expand.grid(
  eta = c(0.1,0.2,0.3),
  max_depth = c(1),
  lambda = c(0.01,0.05,0.1,0.2),
  alpha = c(0.01,0.05,0.1,0.2)
)

results <- data.frame()
# Run grid search
for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  
  model <- xgboost(
    data = dtrain,
    booster = "gbtree",
    objective = "reg:squarederror",
    nrounds = 500,
    eta = params$eta,
    max_depth = params$max_depth,
    lambda = params$lambda,
    alpha = params$alpha,
    #  subsample = 0.7,
    #  colsample_bytree = 0.7,
    verbose = 0
  )
  
  preds <- predict(model, newdata = dtest)
  R2 <- 1 - sum((y_true - preds)^2) / sum((y_true - mean(y_true))^2)
  
  results <- rbind(results, cbind(params, R2))
  print(cbind(params, R2))
}

results <- results[order(-results$R2), ]
print(head(results))     
