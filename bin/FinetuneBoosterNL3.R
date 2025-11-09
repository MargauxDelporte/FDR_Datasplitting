rm(list = ls())
mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
setwd(mywd)

source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Functions/TriangleBoosterTrainMS.R'))

num_split <- 3
p <- 100
p0 <- 10
q <- 0.1
delta <- 10
maximize_power=function(n,myeta,mymax_depth,mylambda,myalpha){
    ### algorithmic settings
  set.seed(456)
  signal_index <- sample(c(1:p), size = p0, replace = F)

  params =list(
    objective = "reg:squarederror",
    eta       = myeta,
    max_depth = mymax_depth,
    lambda    = mylambda,
    alpha     = myalpha
  )
    # simulate data
    n1 <- floor(n/2); n2 <- n - n1
    X1 <- matrix(rnorm(n1*p, mean= 1), n1, p)
    X2 <- matrix(rnorm(n2*p, mean=-1), n2, p)
    X  <- rbind(X1, X2)
    beta_star <- numeric(p)
    beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))*1000
    y <- (X^2 %*% beta_star + rnorm(n))
    
    # run your custom methods
    g1 <- ApplyTriangleBoostTrain( X = X, y = y, q = q, num_split = num_split,param=params,
                                   signal_index = signal_index, myseed = 1)

    # bind all rows
    ResultsDataFrame <- rbind(
      c("Boost DS",n,myeta,mymax_depth,mylambda,myalpha, g1$DS_fdp,g1$DS_power),
      c("Boost MS",n,myeta,mymax_depth,mylambda,myalpha, g1$MDS_fdp,g1$MDS_power)
    )
    return(ResultsDataFrame)
  }
  
param_grid <- expand.grid(
  n=seq(from=500,to=1000,by=50),
  myeta               = c(0.005, 0.01, 0.05, 0.1),
  mymax_depth         = c(6, 8, 10, 12, 15),
  mylambda            = c(0, 0.1, 1),
  myalpha             = c(0, 0.1, 1),
  stringsAsFactors  = FALSE
)

# set up parallel backend
library(doParallel)
library(foreach)

# Set number of cores
num_cores <- 20#parallel::detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Parallelized evaluation i=1
results <- foreach(i = 1:nrow(param_grid), .combine = rbind, .packages = c("xgboost")) %dopar% {
  row <- param_grid[i, ]
  
  maximize_power(
    n            = row$n,
    myeta        = row$myeta,
    mymax_depth  = row$mymax_depth,
    mylambda     = row$mylambda,
    myalpha      = row$myalpha
  )
}

# Stop the cluster
stopCluster(cl)
