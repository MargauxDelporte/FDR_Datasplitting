library(caret)
library(ranger)
library(xgboost)
library(kernlab)
library(earth)
library(caretEnsemble)
#setwd("~/Desktop/temp") j=1 model=lm

permR2EnsembleTrain<-function(data,j,model){
  dataPerm<-data[,-1]
  dataPerm[,j]<-sample(data[,j+1],replace=FALSE)
  predictLM<-predict(model,newdata=as.matrix(dataPerm))
  rsquared=1-sum((data$y-predictLM)^2)/sum((data$y-mean(data$y))^2)
  return(rsquared)
}

myEnsemble=function(){
  rf_fit <- ranger(
    formula = as.formula(paste("y ~", paste(predictors, collapse = "+"))),
    data = df_train2,
    num.trees = 800,
    mtry = max(1, floor(sqrt(p))),
    min.node.size = 5,
    splitrule = "variance",
    sample.fraction = 0.632,
    seed = 456
  )
  
  ## 2b) ExtraTrees (Extremely Randomized Trees)
  et_fit <- ranger(
    formula = as.formula(paste("y ~", paste(predictors, collapse = "+"))),
    data = df_train2,
    num.trees = 800,
    mtry = max(1, floor(p/3)),
    min.node.size = 5,
    splitrule = "extratrees",
    sample.fraction = 1.0,   # often used with ExtraTrees
    seed = 456
  )
  
  ## 2c) XGBoost (gbtree)
  # XGBoost needs a numeric matrix:
  x_train2 <- as.matrix(df_train2[, predictors, drop = FALSE])
  y_train2 <- df_train2$y
  x_val    <- as.matrix(df_val[, predictors, drop = FALSE])
  y_val    <- df_val$y
  x_test   <- as.matrix(df_test[, predictors, drop = FALSE])
  y_test   <- df_test$y
  
  dtrain2 <- xgb.DMatrix(x_train2, label = y_train2)
  dval    <- xgb.DMatrix(x_val,    label = y_val)
  dtest   <- xgb.DMatrix(x_test,   label = y_test)
  
  xgb_params <- list(
    booster = "gbtree",
    objective = "reg:squarederror",
    eval_metric = "rmse",
    eta = 0.05,
    max_depth = 6,
    min_child_weight = 1,
    subsample = 0.8,
    colsample_bytree = 0.8
  )
  
  # Early-stopping using only TRAIN2->VAL (to avoid touching TEST)
  watch <- list(train = dtrain2, val = dval)
  xgb_fit <- xgb.train(
    params = xgb_params,
    data = dtrain2,
    nrounds = 5000,
    watchlist = watch,
    early_stopping_rounds = 100,
    verbose = 0
  )
  
}
ApplyEnsembleTrain<-function(X, y, q,best_nrounds=1000,amountTrain=0.5,amountTest=1-amountTrain,myseed,mybooster='gblinear',num_split=50,signal_index=signal_index){
  set.seed(myseed)
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
    rf   = caretModelSpec(method = "ranger", tuneLength = 5,
                          trControl = ctrl, importance = "impurity"),
    xgb  = caretModelSpec(method = "xgbTree", tuneLength = 8,
                          trControl = ctrl),
    svm  = caretModelSpec(method = "svmRadial", tuneLength = 5,
                          trControl = ctrl),
    knn  = caretModelSpec(method = "knn", tuneLength = 10,
                          trControl = ctrl),
    mars = caretModelSpec(method = "earth", tuneLength = 5,
                          trControl = ctrl)
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
  
  remaining_percent=1-amountTrain
  overlap=max(c(0,amountTest-remaining_percent))
  remaining_index<-c(setdiff(c(1:n),train_index),sample(train_index,size=overlap*n))
  sample_index1 <- sample(x = remaining_index, size = amountTest/2 * n, replace = F)
  sample_index2 <- setdiff(remaining_index, sample_index1)

  predictLM1<-predict(lm,newdata=as.matrix(X[sample_index1,]))
  R2orig1<-1-sum((y[sample_index1]-predictLM1)^2)/sum((y[sample_index1]-mean(y[sample_index1]))^2)
  predictLM2<-predict(lm,newdata=as.matrix(X[sample_index2,]))
  R2orig2<-1-sum((y[sample_index2]-predictLM2)^2)/sum((y[sample_index2]-mean(y[sample_index2]))^2)
  beta1<-sapply(1:ncol(X),function(j) permR2BoostTrain(data[sample_index1,],j,lm))-R2orig1
  beta2<-sapply(1:ncol(X),function(j) permR2BoostTrain(data[sample_index2,],j,lm))-R2orig2
  mirror<-sign(beta1*beta2)*(abs(beta1)+abs(beta2))
  selected_index<-SelectFeatures(mirror,abs(mirror),0.1)
  
  ### number of selected variables
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



