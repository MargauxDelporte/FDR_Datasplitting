### High dimension linear model
rm(list = ls())

mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
#mywd='C:/Users/mde4023/Documents/GitHub/FDR_Datasplitting'

setwd(mywd)

source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Functions/TriangleBoosterTrainMS.R'))
source(paste0(mywd,'/Functions/ApplyGBMKnockoff.R'))
#source(paste0(mywd,'/Functions/TriangleGBMTrainMS.R'))

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
library(neuralnet)

library(glmnet)
library(knockoff)
library(mvtnorm)
library(hdi)

### algorithmic settings
num_split <- 50
n <-750
p <- 100
p0 <- 10
q <- 0.1

#500 10 de OG
#600 1000 teveel
#550 1000 geen verschil
#set.seed(124)(123) i=5
set.seed(456)
signal_index <- sample(c(1:p), size = p0, replace = F)

#eta max_depth lambda alpha   mean_R2
#"Boost MS" "750"  "0.005" "6"  "0"  "0"  "0.181818181818182" "0.9"

###choose the parameters
params =list(
  objective = "reg:squarederror",
  eta       = 0.005,
  max_depth = 6,
  lambda    = 0,
  alpha     = 0
)

#######set up the method for the comparison############# i=10 s=10 num_split=1
Compare_SignalStrength <- function(i, s) {
  set.seed(s)
  delta <- i
  signal_index <- sample(c(1:p), size = p0, replace = F)
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
  # FDR methods
  DS_result      <- DS(          X = X, y = y, q = q, num_split = num_split)
  knockoff_result<- ApplyGBMKnockoff(    X = X, y = y, q = q,param=params)
  BH_result      <- MBHq(        X = X, y = y, q = q, num_split = num_split)
  
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
    data.frame(Method = "Boost DS",                Delta = i, FDP = g1$DS_fdp,    Power = g1$DS_power),
    data.frame(Method = "Boost MS",                Delta = i, FDP = g1$MDS_fdp,   Power = g1$MDS_power),
   data.frame(Method = "DataSplitting",           Delta = i, FDP = DS_result$DS_fdp,  Power = DS_result$DS_power),
    data.frame(Method = "MultipleDataSplitting",   Delta = i, FDP = DS_result$MDS_fdp, Power = DS_result$MDS_power),
    data.frame(Method = "Knockoff",                Delta = i, FDP = knockoff_result$fdp, Power = knockoff_result$power),
    data.frame(Method = "Benjamini–Hochberg (BH)", Delta = i, FDP = BH_result$fdp,     Power = BH_result$power)
  )
  
  return(ResultsDataFrame)
}

Compare_SignalStrength(10,8)

#Compare_SignalStrength(7,7)
#######run the code#############
#Results=data.frame()
#for(s in 1:25){
#  for(i in seq(from=5,to=13,by=1)){
#  Results=rbind(Results,Compare_SignalStrength(i,s))
#  print(s)
#  }
#  print(Results)
#  }

library(parallel)

mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
setwd(mywd)

# Source helper and method files

source(file.path(mywd, 'Functions', 'TriangleBoosterTrainMS.R'))
source(file.path(mywd, 'Functions', 'HelperFunctions.R'))
source(file.path(mywd, 'Functions', 'ApplyGBMKnockoff.R'))

# Dai’s routines
#source(file.path(mywd, 'Functions Dai', 'knockoff.R'))
source(file.path(mywd, 'Functions Dai', 'analysis.R'))
source(file.path(mywd, 'Functions Dai', 'MBHq.R'))
source(file.path(mywd, 'Functions Dai', 'DS.R'))
source(file.path(mywd, 'Functions Dai', 'fdp_power.R'))

# Load required packages
pkgs <- c('xgboost','gbm','ranger','MASS','glmnet','knockoff','mvtnorm','hdi',
          'foreach','doParallel')
lapply(pkgs, library, character.only = TRUE)

# === PARAMETER GRID ===
param_grid <- expand.grid(
  s = 1:50,
  i = seq(from = 7, to = 13, by = 1)
)

# === SET UP PARALLEL BACKEND ===
cl <- makeCluster(20)
# export working dir so workers can source
clusterExport(cl, 'mywd')
# have each worker source & load libraries
clusterEvalQ(cl, {
  setwd(mywd)
  source(file.path(mywd, 'Functions', 'TriangleBoosterTrainMS.R'))
  source(file.path(mywd, 'Functions', 'ApplyGBMKnockoff.R'))
  source(file.path(mywd, 'Functions', 'HelperFunctions.R'))
  source(file.path(mywd, 'Functions Dai', 'knockoff.R'))
  source(file.path(mywd, 'Functions Dai', 'analysis.R'))
  source(file.path(mywd, 'Functions Dai', 'MBHq.R'))
  source(file.path(mywd, 'Functions Dai', 'DS.R'))
  source(file.path(mywd, 'Functions Dai', 'fdp_power.R'))
  lapply(c('xgboost','gbm','ranger','MASS','glmnet','knockoff','mvtnorm','hdi'),
         library, character.only = TRUE)
})
registerDoParallel(cl)

# === RUN IN PARALLEL AND WRITE OUT ===
results_list <- foreach(
  k = seq_len(nrow(param_grid)),
  .packages = pkgs,
  .combine  = rbind
) %dopar% {
  s_val <- param_grid$s[k]
  i_val <- param_grid$i[k]
  
  # compute chunk of results
  chunk <- Compare_SignalStrength(i = i_val, s = s_val)
  
  # write out this chunk immediately
  fname <- sprintf("Results_s%02d_i%02d.csv", s_val, i_val)
  write.csv(chunk, file = paste0(mywd,"/Temp/",fname), row.names = FALSE)
  
  # return for final binding
  chunk
}

# === CLEANUP AND FINAL SAVE ===
stopCluster(cl)
warnings()
# combine all and save full dataset
Results <- results_list
write.csv(Results, file = "ResultsNonlinearScenario2.csv", row.names = FALSE)


# === Add 25 additional simulations ===

# Path to your folder
csv_dir <- "C:/Users/mde4023/Downloads/FDR_Datasplitting/Temp"
csv_files <- list.files(
  path       = csv_dir,
  pattern    = "\\.csv$",
  full.names = TRUE
)
# Get full paths of all .csv files
# Read each file into a list of data.frames
data_list <- lapply(csv_files, read.csv, stringsAsFactors = FALSE)
library(dplyr)
all_data <- bind_rows(data_list, .id = "source_file")

##########visualise the results###########
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readr)
mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
mywd <- paste0(mywd,'/Results')

colors <- c("#000000","#FF00FF","#009900", "#99ccff", "#0000FF", "#FF0000")
Results2=Results
names(Results2)=c('Method','SignalStrength','FDR','Power')
Results2$FDR=round(as.numeric(Results2$FDR),3)
Results2$Power=round(as.numeric(Results2$Power),2)
resultsagg <- Results2 %>%
  group_by(Method, SignalStrength) %>%
  summarize(
    Avg_FDR = mean(FDR),
    Avg_Power = mean(Power)
  )
resultsagg$Signal_noisy <- as.numeric(resultsagg$SignalStrength) + runif(nrow(resultsagg), -0.2, 0.2)
View(resultsagg)
resultsagg <- resultsagg %>%
  mutate(Method  = case_when(
    Method == "BH" ~ "Benjamini–Hochberg",
    Method == "DataSplitting" ~ "Dai (single split)",
    Method == "MultipleDataSplitting" ~ "Dai (50 splits)",
    Method == "Knockoff" ~ "Knockoff",
    Method == "Boost DS" ~ "Delporte (single split)",
    Method == "Boost MS" ~ "Delporte (50 splits)",
    TRUE ~ Method  # default if none match
  ))
PowerPlot <- ggplot(resultsagg, aes(x = Signal_noisy, y = as.numeric(Avg_Power), color = Method)) +
  geom_point(size = 3) +
  geom_line()+
  labs(x = "Signal", y = "Power") +
  scale_x_continuous(breaks = seq(from = 5, to = 13, by = 1)) +
  geom_hline(yintercept = 0.8) +
  scale_color_manual(values = colors)+
  coord_cartesian(ylim = c(-0.01, 1)) 
FDRPlot=ggplot(resultsagg, aes(x = Signal_noisy, y = as.numeric(Avg_FDR ), color = Method)) +
  geom_point(size = 3) +geom_line()+
  labs(x = "Signal", y = "FDR")+
  scale_x_continuous(breaks=seq(from=5,to=13,by=1))+
  geom_hline(yintercept=0.1)+
  scale_color_manual(values = colors)
PlotPermute=ggarrange(
  PowerPlot, FDRPlot,
  common.legend = TRUE, legend = "right"
)
PlotPermute
ggsave("NLScenario.png",
       plot   = PlotPermute,
       width  = 8,
       height = 8/18*8,
       units  = "in",
       dpi    = 100)
