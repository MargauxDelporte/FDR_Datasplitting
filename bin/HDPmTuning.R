library(caret)


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
y <- scale(X %*% beta_star + rnorm(n))



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

# crossvalidation lasso
cv_lasso <- cv.glmnet(x = as.matrix(Xtrain), y = y[train_index], alpha = 1)
        
best_lambda <- cv_lasso$lambda.min
        
# Fit final model with best lambda
model <- glmnet(x = as.matrix(Xtrain), y = y[train_index], alpha = 1, lambda = best_lambda)
        
# Predict on test set
preds <- predict(model, newx = data.matrix(X[sample_index1, ]))
        
# Compute R² on test set
y_true <- y[sample_index1]
R2 <- 1 - sum((y_true - preds)^2) / sum((y_true - mean(y_true))^2)
  
ApplyTriangleLassoHD(X = X, y = y, q = q, num_split = num_split,signal_index = signal_index, myseed = 1)     
