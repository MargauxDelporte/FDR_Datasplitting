### High dimension linear model
rm(list = ls())

mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
#mywd='C:/Users/mde4023/Documents/GitHub/FDR_Datasplitting'
#mywd='C:/Users/marga/Downloads/FDR_Datasplitting'
setwd(mywd)

source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Functions/TriangleBoosterTrainMS.R'))
source(paste0(mywd,'/Functions/ApplyGBMKnockoff.R'))
source('C:/Users/mde4023/Downloads/FDR_Datasplitting/Results/ResultsHDNL Scenario/MarsParallelHD.R')

source(paste0(mywd,'/Functions Dai/knockoff.R'))
source(paste0(mywd,'/Functions Dai/analysis.R'))
source(paste0(mywd,'/Functions Dai/MBHq.R'))
source(paste0(mywd,'/Functions Dai/DS.R'))
source(paste0(mywd,'/Functions Dai/fdp_power.R'))

#devtools::install_github("Jeremy690/DSfdr/DSfdr",force = TRUE)
library(xgboost)
library(gbm)
library(ranger)
library(MASS)

library(glmnet)
library(knockoff)
library(mvtnorm)
library(hdi)
#C:/Users/mde4023/Downloads/FDR_Datasplitting/Results/ResultsHDNL Scenario/HDNonlinearScenario6_mars.R
#algorithmic settings
num_split <- 50
n <-500
p <- 550
p0 <- 10#25
q <- 0.10
amountTest=0.5
amountTrain=0.5
set.seed(11212025)
signal_index <- sample(c(1:p), size = p0, replace = F)
###choose the parameters
params =list(
  objective = "reg:squarederror",
  eta       = 0.005,
  max_depth = 6,
  lambda    = 0,
  alpha     = 0
)

#######set up the method for the comparison############# i=7;s=1 num_split=1
Compare_SignalStrength <- function(i, s,other=T) {
  set.seed(s)
  delta <- i
  # simulate data
  n1 <- floor(n/2); n2 <- n - n1
  X1 <- matrix(rnorm(n1*p, mean= 1), n1, p)
  X2 <- matrix(rnorm(n2*p, mean=-1), n2, p)
  X  <- rbind(X1, X2)
  beta_star <- numeric(p)
  beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))*1000
  y <- (X^2 %*% beta_star+ rnorm(n))
  # run your custom methods
  g1 <- ApplyMarsTrain_HDparallel( X = X, y = y, q = q, num_split = num_split,signal_index = signal_index, myseed = 1)
  # FDR methods
  if(other){
    DS_result      <- DS(          X = X, y = y, q = q, num_split = num_split)
    knockoff_result<- ApplyGBMKnockoff(    X = X, y = y, q = q,param=params)
    BH_result      <- MBHq(        X = X, y = y, q = q, num_split = num_split)
  }
  # init empty results df
  ResultsDataFrame <- data.frame(
    Method = character(),
    Delta  = numeric(),
    FDP    = numeric(),
    Power  = numeric(),
    stringsAsFactors = FALSE
  )
  
  # bind all rows
  ResultsDataFrame <- rbind(
    ResultsDataFrame,
    data.frame(Method = "Mars DS",                Delta = i, FDP = g1$DS_fdp,    Power = g1$DS_power),
    data.frame(Method = "Mars MS",                Delta = i, FDP = g1$MDS_fdp,   Power = g1$MDS_power))
  if(other){
    ResultsDataFrame <- rbind(
      ResultsDataFrame,
      data.frame(Method = "DataSplitting",           Delta = i, FDP = DS_result$DS_fdp,  Power = DS_result$DS_power),
      data.frame(Method = "MultipleDataSplitting",   Delta = i, FDP = DS_result$MDS_fdp, Power = DS_result$MDS_power),
      data.frame(Method = "Knockoff",                Delta = i, FDP = knockoff_result$fdp, Power = knockoff_result$power),
      data.frame(Method = "Benjamini–Hochberg (BH)", Delta = i, FDP = BH_result$fdp,     Power = BH_result$power)
    )
  }
  return(ResultsDataFrame)
}
Compare_SignalStrength(7,1,F)

#Compare_SignalStrength(13,13,F)
# build grid
param_grid <- expand.grid(
  s = 1:50,
  i = 8:13
)

# make sure output dir exists
mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting/Results/ResultsHDNL Scenario/'
out_dir <- file.path(mywd, "Temp2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# iterate over ROWS, not columns
for (k in seq_len(nrow(param_grid))) {
  s_val <- param_grid$s[k]
  i_val <- param_grid$i[k]
  
  # compute chunk (wrap in tryCatch to avoid breaking the whole run)
  chunk <- tryCatch(
    Compare_SignalStrength(i = i_val, s = s_val),
    error = function(e) {
      message(sprintf("Failed at s=%s, i=%s: %s", s_val, i_val, conditionMessage(e)))
      return(NULL)
    }
  )
  if (is.null(chunk)) next
  if (!is.data.frame(chunk)) chunk <- as.data.frame(chunk)
  
  # write out this chunk immediately
  fname <- sprintf("Results_s%02d_i%02d.csv", s_val, i_val)
  write.csv(chunk, file = file.path(out_dir, fname), row.names = FALSE)
}

# === CLEANUP AND FINAL SAVE ===

# Path to your folder
# === CLEANUP AND FINAL SAVE ===

# Path to your folder
csv_dir <- "C:/Users/mde4023/Downloads/FDR_Datasplitting/Temp2"
csv_files <- list.files(
  path       = csv_dir,
  pattern    = "\\.csv$",
  full.names = TRUE
)
warnings()


