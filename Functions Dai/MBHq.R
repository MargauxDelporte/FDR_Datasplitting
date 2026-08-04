library(glmnet)

linearModel <- function(x, y,
                            family = "gaussian",
                            nfolds = 10,
                            s = "lambda.1se",
                            standardize = TRUE) {
  
  x <- as.matrix(x)
  y <- as.numeric(y)
  
  fit <- lm(y~x)
  pvals=summary(fit)$coef[-1,4]
  return(pvals)
}

MBHq <- function(X, y, q){
  pvalues = linearModel(X, y)
  result=p.adjust(pvalues, method = "BY")
  selected_index=which(result<q)
  result=CalculateFDP_Power(selected_index, signal_index)
  BH_fdp <- result$fdp
  BH_power <- result$power
  return(list(fdp = BH_fdp, power = BH_power))
}
