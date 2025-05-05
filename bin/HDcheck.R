library(MASS)
library(xgboost)
mywd='C:/Users/mde4023/Documents/GitHub/FDR_Datasplitting'

setwd(mywd)

source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Functions/TriangleBoosterTrainMS.R'))
num_split <- 10
n <-1500
p <- 2000
p0 <- 50
q <- 0.1
delta <-5
signal_index <- sample(c(1:p), size = p0, replace = F)
# simulate data
X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
beta_star <- numeric(p)
beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))*100
y <- scale(X %*% beta_star + rnorm(n))

# run your custom methods
g1 <- ApplyTriangleBoostTrain( X = X, y = y, q = q, num_split = num_split,mybooster='gblinear',
                               signal_index = signal_index, myseed = 1)
g1
