### High dimension linear model
rm(list = ls())
getwd()
#mywd='C:/Users/mde4023/OneDrive - Weill Cornell Medicine/0 Projects/FDR_Datasplitting'
mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
#mywd='~/FDR_Datasplitting'

setwd(mywd)

source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Functions/TriangleLassoHD.R'))


source(paste0(mywd,'/Functions Dai/knockoff.R'))
source(paste0(mywd,'/Functions Dai/analysis.R'))
source(paste0(mywd,'/Functions Dai/MBHq.R'))
source(paste0(mywd,'/Functions Dai/DS.R'))
source(paste0(mywd,'/Functions Dai/fdp_power.R'))

#devtools::install_github("Jeremy690/DSfdr/DSfdr",force = TRUE)
library(glmnet)
library(knockoff)
library(mvtnorm)
library(hdi)
library(MASS)
### algorithmic settings
num_split <- 10
n <-1500
p <- 2000
p0 <- 25
q <- 0.1

#set.seed(124)(123) i=5
set.seed(456)
signal_index <- sample(c(1:p), size = p0, replace = F)
#######set up the method for the comparison############# i=10
Compare_SignalStrength <- function(i, s) {
  set.seed(s)
  delta <- i
  
  # simulate data
  X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
  beta_star <- numeric(p)
  beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))
  y <- scale(X %*% beta_star + rnorm(n))
  
  # run your custom methods
  g1 <- ApplyTriangleLassoHD(X = X, y = y, q = q, num_split = num_split,
                             signal_index = signal_index, myseed = 1)     
  
  # FDR methods
  DS_result      <- DS(          X = X, y = y, q = q, num_split = num_split)
  knockoff_result<- knockoff(    X = X, y = y, q = q)
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
    data.frame(seed=s, Method = "Boost DS",                Delta = i, FDP = g1$DS_fdp,    Power = g1$DS_power),
    data.frame(seed=s, Method = "Boost MS",                Delta = i, FDP = g1$MDS_fdp,   Power = g1$MDS_power),
    data.frame(seed=s, Method = "DataSplitting",           Delta = i, FDP = DS_result$DS_fdp,  Power = DS_result$DS_power),
    data.frame(seed=s, Method = "MultipleDataSplitting",   Delta = i, FDP = DS_result$MDS_fdp, Power = DS_result$MDS_power),
    data.frame(seed=s, Method = "Knockoff",                Delta = i, FDP = knockoff_result$fdp, Power = knockoff_result$power),
    data.frame(seed=s, Method = "Benjamini–Hochberg (BH)", Delta = i, FDP = BH_result$fdp,     Power = BH_result$power)
  )
  print(ResultsDataFrame)
  return(ResultsDataFrame)
}


library(parallel)

# Source helper and method files

source(file.path(mywd, 'Functions', 'TriangleLassoHD.R'))
source(file.path(mywd, 'Functions', 'HelperFunctions.R'))

# Dai’s routines
source(file.path(mywd, 'Functions Dai', 'knockoff.R'))
source(file.path(mywd, 'Functions Dai', 'analysis.R'))
source(file.path(mywd, 'Functions Dai', 'MBHq.R'))
source(file.path(mywd, 'Functions Dai', 'DS.R'))
source(file.path(mywd, 'Functions Dai', 'fdp_power.R'))

# Load required packages
pkgs <- c('MASS','glmnet','knockoff','mvtnorm','hdi',
          'foreach','doParallel')
lapply(pkgs, library, character.only = TRUE)

# === PARAMETER GRID ===
param_grid <- expand.grid(
  s = c(1:6,26:50),
  i = seq(from = 7, to = 13, by = 1)
)
param_grid

# === SET UP PARALLEL BACKEND ===
cl <- makeCluster(15)#detectCores() - 4)
# export working dir so workers can source
clusterExport(cl, 'mywd')
# have each worker source & load libraries
clusterEvalQ(cl, {
  setwd(mywd)
  source(file.path(mywd, 'Functions', 'TriangleLassoHD.R'))
  source(file.path(mywd, 'Functions', 'HelperFunctions.R'))
  source(file.path(mywd, 'Functions Dai', 'knockoff.R'))
  source(file.path(mywd, 'Functions Dai', 'analysis.R'))
  source(file.path(mywd, 'Functions Dai', 'MBHq.R'))
  source(file.path(mywd, 'Functions Dai', 'DS.R'))
  source(file.path(mywd, 'Functions Dai', 'fdp_power.R'))
  lapply(c('MASS','glmnet','knockoff','mvtnorm','hdi','foreach','doParallel'),
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


param_grid <- expand.grid(
  s = c(1:6,
  i = seq(from = 7, to = 13, by = 1)
)
k=1
results=c()
for(k in 1:nrow(param_grid)){
  # compute chunk of results
  s_val <- param_grid$s[k]
  i_val <- param_grid$i[k]
  fname <- sprintf("C:/Users/mde4023/Downloads/FDR_Datasplitting/Results HD/Results_s%02d_i%02d.csv", s_val, i_val)
  rf=read.csv(file = fname)
  results=rbind(results,rf)
}
library(openxlsx)
#write.xlsx(results,'C:/Users/mde4023/Downloads/FDR_Datasplitting/Results/results_HD.xlsx',rowNames=F)

##########visualise the results###########
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readr)

Results=results
#Results=read.xlsx('VanillaResults.xlsx')
colors <- c("#000000","#FF00FF","#009900", "#99ccff", "#0000FF", "#FF0000")
Results2=Results
names(Results2)=c('seed','Method','SignalStrength','FDP','Power')
Results2$FDP=round(as.numeric(Results2$FDP),3)
Results2$Power=round(as.numeric(Results2$Power),2)
resultsagg <- Results2 %>%
  group_by(Method, SignalStrength) %>%
  summarize(
    Avg_FDR = mean(FDP),
    Avg_Power = mean(Power)
  )
resultsagg$Signal_noisy <- as.numeric(resultsagg$SignalStrength) + runif(nrow(resultsagg), -0.2, 0.2)

resultsagg <- resultsagg %>%
  mutate(Method  = case_when(
    Method == "BH" ~ "Benjamini–Hochberg",
    Method == "DataSplitting" ~ "Dai (single split)",
    Method == "MultipleDataSplitting" ~ "Dai (10 splits)",
    Method == "Knockoff" ~ "Knockoff",
    Method == "Boost DS" ~ "Delporte (single split)",
    Method == "Boost MS" ~ "Delporte (10 splits)",
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
  geom_point(size = 3) +
  geom_line()+
  labs(x = "Signal", y = "FDP")+
  scale_x_continuous(breaks=seq(from=5,to=13,by=1))+
  geom_hline(yintercept=0.1)+
  scale_color_manual(values = colors)
PlotPermute=ggarrange(
  PowerPlot, FDRPlot,
  common.legend = TRUE, legend = "right"
)

ggsave("HDScenario.png",
       plot   = PlotPermute,
       width  = 8,
       height = 8/18*8,
       units  = "in",
       dpi    = 100)
