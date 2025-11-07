library(caret)
library(ranger)
library(xgboost)
library(kernlab)
library(earth)
library(caretEnsemble)
library(doParallel)


#setwd("~/Desktop/temp") j=1 model=lm

permR2EnsembleTrain<-function(data,j,model){
  dataPerm<-data[,-1]
  dataPerm[,j]<-sample(data[,j+1],replace=FALSE)
  predictLM<-predict(model,newdata=as.matrix(dataPerm))
  rsquared=1-sum((data$y-predictLM)^2)/sum((data$y-mean(data$y))^2)
  return(rsquared)
}

## 3) Tight, smart tune grids (instead of big tuneLength sweeps)

# ranger: fast, strong baseline
grid_rf <- expand.grid(
  mtry       = pmax(1, round(c(0.05, 0.15, 0.3) * (ncol(dataTrain) - 1))),
  splitrule  = "gini",
  min.node.size = c(1, 5, 10)
)

# xgbTree: keep nrounds modest; subsample/colsample speed up a lot
grid_xgb <- expand.grid(
  nrounds = c(100, 200),
  max_depth = c(3, 5),
  eta = c(0.05, 0.1),
  gamma = 0,
  colsample_bytree = c(0.6, 0.8),
  min_child_weight = c(1, 3),
  subsample = c(0.7, 0.9)
)
# also make sure XGBoost uses threads
set.seed(1)
xgb.set.config(verbosity = 0) # optional
# caret passes threads via env var; also set:
Sys.setenv(OMP_NUM_THREADS = cores)

# SVM: small grid; cache size helps; scale beforehand
grid_svm <- expand.grid(
  C     = c(0.5, 1, 2),
  sigma = c(0.01, 0.05, 0.1)
)

# KNN: keep k small set; gets slow with big p, so consider dropping if p is large
grid_knn <- data.frame(k = c(3, 5, 9, 15))

# MARS: limit basis functions (nk) to speed up
grid_mars <- expand.grid(
  degree = c(1, 2),
  nprune = c(10, 20, 35)
)

ApplyEnsembleTrain<-function(X, y, q,best_nrounds=1000,amountTrain=0.5,amountTest=1-amountTrain,myseed,mybooster='gblinear',num_split=50,signal_index=signal_index){
  set.seed(myseed)
  data<-data.frame(cbind(y,X))
  n <- dim(X)[1]; p <- dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  fdp <- rep(0, num_split)
  power <- rep(0, num_split)
  num_select <- rep(0, num_split)
  data<-data.frame(cbind(y,X))
  #for(iter in 1:num_split){  }
  train_index<-sample(x = c(1:n), size = amountTrain * n, replace = F)
  dataTrain<-data[train_index,]
  colnames(dataTrain)<-c('y',paste0('X',1:p))
  colnames(data)<-c('y',paste0('X',1:p))
  
  # ----------------------------
  # 2. Define resampling + tuning controls
  # ----------------------------
  ctrl <- trainControl(
    method = "cv", number = 5,
    savePredictions = "final",
    allowParallel = TRUE
  )
  
  # ----------------------------
  # 3. Specify a suite of nonlinear learners
  # ----------------------------
  model_list <- list(
    rf   = caretModelSpec(
      method = "ranger",
      tuneGrid = grid_rf,
      trControl = ctrl,
      importance = "none",         # turn off to save time
      num.trees = 500,             # 500 is often enough; 1000+ slows down
      respect.unordered.factors = "order",
      num.threads = cores
    ),
    xgb  = caretModelSpec(
      method = "xgbTree",
      tuneGrid = grid_xgb,
      trControl = ctrl,
      verbose = FALSE              # quiet logs
    ),
    mars = caretModelSpec(
      method = "earth",
      tuneGrid = grid_mars,
      trControl = ctrl,
      pmethod = "none",            # no internal CV in earth (faster)
      nk = 25                      # cap basis funcs
    )
  )
  
  # ----------------------------
  # 4. Train all models with caretList (needs caretEnsemble)
  # ----------------------------
  models <- caretList(
    y ~ ., data = dataTrain,
    trControl = ctrl,
    tuneList = model_list,
    continue_on_fail = TRUE
  )
  # ----------------------------
  # 5. Stack them with a meta-learner
  # ----------------------------
  stack_ctrl <- trainControl(method = "cv", number = 5, savePredictions = "final")
  stack_model <- caretStack(
    models,
    method = "glmnet",           # regularized linear blender
    metric = "RMSE",
    trControl = stack_ctrl
  )
  
  # ----------------------------
  # 6. Use the stacked model like your previous lm
  # ----------------------------
  lm <- stack_model
  
  remaining_index<-c(setdiff(c(1:n),train_index))
  sample_index1 <- sample(x = remaining_index, size = amountTest/2 * n, replace = F)
  sample_index2 <- setdiff(remaining_index, sample_index1)

  predict_TRAIN<-predict(lm,newx=as.matrix(X[train_index,]))
  R2orig_TRAIN<-1-sum((y[train_index]-predict_TRAIN)^2)/sum((y[train_index]-mean(y[train_index]))^2)
  R2orig_TRAIN
  
  predictLM1<-predict(lm,newx=as.matrix(X[sample_index1,]))
  predictLM2<-predict(lm,newx=as.matrix(X[sample_index2,]))
  
  R2orig1<-1-sum((y[sample_index1]-predictLM1)^2)/sum((y[sample_index1]-mean(y[sample_index1]))^2)
  R2orig2<-1-sum((y[sample_index2]-predictLM2)^2)/sum((y[sample_index2]-mean(y[sample_index2]))^2)
  
  result=c(R2orig_TRAIN,R2orig1,R2orig2)
  result

  return(result)
}



