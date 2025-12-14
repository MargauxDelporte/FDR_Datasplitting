### High dimension linear model
rm(list = ls())

mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
#mywd='C:/Users/mde4023/Documents/GitHub/FDR_Datasplitting'
#C:\Users\mde4023\Downloads\FDR_Datasplitting

setwd(mywd)
source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Functions/TriangleLinRegTrainMS.R'))
source(paste0(mywd,'/Functions Dai/knockoff.R'))
source(paste0(mywd,'/Functions Dai/analysis.R'))
source(paste0(mywd,'/Functions Dai/MBHq.R'))
source(paste0(mywd,'/Functions Dai/DS.R'))
source(paste0(mywd,'/Functions Dai/fdp_power.R'))

#devtools::install_github("Jeremy690/DSfdr/DSfdr",force = TRUE)

library(MASS)
library(glmnet)
library(knockoff)
library(mvtnorm)
library(hdi)
library(dplyr)
library(ggplot2)
library(ggpubr)


### algorithmic settings
num_split <- 50
n <-1500
p <- 250
p0 <- 25
q <- 0.1
#set.seed(124)(123) i=5
set.seed(456)
signal_index <- sample(c(1:p), size = p0, replace = F)

################################Triange: train test test each 30####################################
Compare_SignalStrength=function(i,s){
  set.seed(s)
  ResultsDataFrame=data.frame()
  delta <- i
  X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
  
  ### randomly generate the true beta i=4
  beta_star <- rep(0, p)
  beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
  
  ### generate y
  y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
  
  ###my own methods:
  g=ApplyTriangleLinRegTrain(X=as.data.frame(X), y, q=0.1,num_split=num_split, signal_index=signal_index, amountTrain=0.333, myseed = 1)
  
  ResultsDataFrame=c('LinReg DS',i, as.numeric(g$DS_fdp),as.numeric(g$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('LinReg MS',i, as.numeric(g$MDS_fdp),as.numeric(g$MDS_power)))
  
  ### Competition
  DS_result <- DS(X,y, num_split=num_split, q=0.1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('DataSplitting',i,DS_result$DS_fdp,DS_result$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('MultipleDataSplitting',i,DS_result$MDS_fdp,DS_result$MDS_power))
  
  knockoff_result <- knockoff(X, y, q=0.1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('Knockoff',i,knockoff_result$fdp,knockoff_result$power))
  
  BH_result <- MBHq(X, y, q=0.1, num_split)
  ResultsDataFrame=rbind(ResultsDataFrame,c('BH',i,BH_result$fdp,BH_result$power))
  
  ### save data
  return(ResultsDataFrame)}

library(parallel)
mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
setwd(mywd)

# Source helper and method files
source(file.path(mywd, 'Functions','HelperFunctions.R'))
source(file.path(mywd, 'Functions', 'TriangleLinRegTrainMS.R'))

# Dai’s routines
source(file.path(mywd, 'Functions Dai', 'knockoff.R'))
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
cl <- makeCluster(5)
# export working dir so workers can source
clusterExport(cl, 'mywd')
# have each worker source & load libraries
clusterEvalQ(cl, {
  setwd(mywd)
  source(file.path(mywd, 'Functions', 'HelperFunctions.R'))
  source(file.path(mywd, 'Functions', 'TriangleLinRegTrainMS.R'))
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
  write.csv(chunk, file = paste0(mywd,"/Results/ResultsLinearScenario/Temp/",fname), row.names = FALSE)
  
  # return for final binding
  chunk
}
# === CLEANUP AND FINAL SAVE ===
stopCluster(cl)
warnings()
# combine all and save full dataset
Results <- results_list
#write.csv(Results, file = "All_Results_vanillaLinear.csv", row.names = FALSE)
Results=as.data.frame(Results)
names(Results)=c('Method','SignalStrength', 'FDR','Power')
head(Results)
library(openxlsx)
library(readr)
write.xlsx(Results,file=paste0(mywd,"/Results/VanillaScenarioLinear.xlsx"))

names(Results)=c('Method','SignalStrength', 'FDR','Power')

####Visualise the results
# Create a color palette based on the number of unique methods
Results2=Results
#Results=read.xlsx('VanillaResults.xlsx')
colors <- c("#FF0000", "#00FF00", "#0000FF", "#000000", "#FF00FF", )
sort(unique(Results2$Method))
colors <- c("#000000","#FF00FF","#009900", "#99ccff", "#0000FF", "#FF0000")
Results2=Results
Results2$FDR=round(as.numeric(Results2$FDR),3)
Results2$Power=round(as.numeric(Results2$Power),2)
#Results2=subset(Results2,SignalStrength%in%as.character(7:13))
resultsagg <- Results2 %>%
  group_by(Method, SignalStrength) %>%
  summarize(
    Avg_FDR = mean(FDR),
    Avg_Power = mean(Power)
  )
View(resultsagg)
resultsagg$Signal_noisy <- as.numeric(resultsagg$SignalStrength) + runif(nrow(resultsagg), -0.2, 0.2)

resultsagg <- resultsagg %>%
  mutate(Method  = case_when(
    Method == "BH" ~ "Benjamini–Hochberg (BH)",
    Method == "DataSplitting" ~ "Dai (single split)",
    Method == "MultipleDataSplitting" ~ "Dai (50 splits)",
    Method == "Knockoff" ~ "Knockoff",
    Method == "LinReg DS" ~ "Delporte (single split)",
    Method == "LinReg MS" ~ "Delporte (50 splits)",
    TRUE ~ Method  # default if none match
  ))
PowerPlot <- ggplot(resultsagg, aes(x = Signal_noisy, y = as.numeric(Avg_Power), color = Method)) +
  geom_point(size = 3) +
  geom_line()+
  labs(x = "Signal", y = "Power") +
  scale_x_continuous(breaks = seq(from = 3, to = 13, by = 1)) +
  geom_hline(yintercept = 0.8) +
  scale_color_manual(values = colors)+
  coord_cartesian(ylim = c(0.1, 1)) 
FDRPlot=ggplot(resultsagg, aes(x = Signal_noisy, y = as.numeric(Avg_FDR ), color = Method)) +
  geom_point(size = 3) +
  geom_line()+
  labs(x = "Signal", y = "FDR")+
  scale_x_continuous(breaks=seq(from=3,to=13,by=1))+
  geom_hline(yintercept=0.1)+
  scale_color_manual(values = colors)
PlotPermute=ggarrange(
  PowerPlot, FDRPlot,
  common.legend = TRUE, legend = "right"
)
PlotPermute

8/18*8
ggsave("VanillaScenario.png",
       plot   = PlotPermute,
       width  = 8,
       height = 8/18*8,
       units  = "in",
       dpi    = 100)
