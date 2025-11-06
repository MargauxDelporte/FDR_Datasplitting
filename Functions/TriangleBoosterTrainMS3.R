params_best <- list(
  objective = "reg:squarederror",
  max_depth = 4,          # allow higher-order interactions
  min_child_weight = 1,
  eta = 0.05,             # smaller step size
  subsample = 0.8,
  colsample_bytree = 0.8,
  reg_lambda = 1.0,
  reg_alpha = 0.0,
  nthread = 0             # use all CPU threads
)
best_nrounds <- 88L

permR2TriangleBoostTrain2<-function(data,Y,j,model){
  Xperm <- data
  Xperm[, j] <- sample(data[, j], replace = FALSE)
  pred_perm <- predict(model, Xperm)
  rsq_perm <- 1 - sum((Y - pred_perm)^2) / sum((Y - mean(Y))^2)
  rsq_perm
    return(rsq_perm)
}

ApplyTriangleBoostTrain2<-function(X, y, q=0.1,myseed=1,num_split=1,signal_index=signal_index){
  set.seed(myseed)
 # amountTrain=0.333
#  amountTest=1-amountTrain
  n_train <- floor(0.6 * n)
  n_valid <- floor(0.2 * n)
  
  n <- dim(X)[1]; p <- dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  fdp <- rep(0, num_split)
  power <- rep(0, num_split)
  num_select <- rep(0, num_split)

  for(iter in 1:num_split){
  idx <- sample.int(n, n)
  id_tr <- idx[1:n_train]
  id_va <- idx[(n_train + 1):(n_train + n_valid)]
  id_te <- idx[(n_train + n_valid + 1):n]
    
  Xtr <- X[id_tr, , drop = FALSE]; ytr <- y[id_tr]
  Xva <- X[id_va, , drop = FALSE]; yva <- y[id_va]
  Xte <- X[id_te, , drop = FALSE]; yte <- y[id_te]
  
  dtr   <- xgb.DMatrix(Xtr, label = ytr)
  dtva <- xgb.DMatrix(Xva, label =  yva)
  dte   <- xgb.DMatrix(Xte, label = yte)
  
  fit_best <- xgb.train(
    params = params_best,
    data = dtr,
    nrounds = 2000,                 # generous cap
    verbose = 0
  )
  
  pred_xgb0 <- predict(fit_best, dtr)
  pred_xgb1 <- predict(fit_best, dtva)
  pred_xgb2 <- predict(fit_best, dte)



  R2orig_TRAIN<-1-sum((ytr-pred_xgb0)^2)/sum((ytr-mean(ytr))^2)
  R2orig_TRAIN
  
  R2orig1<-1-sum((yva-pred_xgb1)^2)/sum((yva-mean(yva))^2)
  R2orig2<-1-sum((yte-pred_xgb2)^2)/sum((yte-mean(yte))^2)
  
  Rnew1<-sapply(1:ncol(X),function(j) permR2TriangleBoostTrain2(data=Xva,Y=yva,j,model=fit_best))
  Rnew2<-sapply(1:ncol(X),function(j) permR2TriangleBoostTrain2(data=Xte,Y=yte,j,fit_best))

  beta1=R2orig1-Rnew1
  beta2=R2orig2-Rnew2

  mirror<-sign(beta1*beta2)*(abs(beta1)+abs(beta2))
  hist(mirror[-signal_index])
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


