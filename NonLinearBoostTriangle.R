### High dimension linear model
rm(list = ls())

mywd='C:/Users/mde4023/OneDrive - Weill Cornell Medicine/0 Projects/FDR_Datasplitting'
mywd='C:/Users/mde4023/Documents/GitHub/FDR_Datasplitting'

setwd(mywd)
source('HelperFunctions.R')
source('TriangleBoosterTrainMS.R')
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
#library(DSfdr)
library(ggplot2)
# library(igraph)
library(RColorBrewer)
library(openxlsx)
library(ggpubr)
library(viridis)

### algorithmic settings
num_split <- 10
n <-1500
p <- 250
p0 <- 25
q <- 0.1
#set.seed(124)(123) i=5
set.seed(456)
signal_index <- sample(c(1:p), size = p0, replace = F)

################################Triange: train test test each 30####################################
Compare_SignalStrenght=function(i,s){
  set.seed(s)
  ResultsDataFrame=data.frame()
  delta <- i
  X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
  
  ### randomly generate the true beta i=40
  beta_star <- rep(0, p)
  beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
  
  ### generate y
  y <- X^2%*%beta_star + rnorm(n, mean = 0, sd = 1)
  
  ###my own methods:
  g=ApplyTriangleBoostTrain(X=as.data.frame(X), y, q=0.1,num_split=10,,mybooster='gbtree',signal_index=signal_index, amountTrain=0.333, myseed = 1)
  
  ResultsDataFrame=c('TreeBoost DS',i, as.numeric(g$DS_fdp),as.numeric(g$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('TreeBoost MS',i, as.numeric(g$MDS_fdp),as.numeric(g$MDS_power)))
  
  ### Competition
  DS_result <- DS(X,y, num_split=10, q=0.1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('DataSplitting',i,DS_result$DS_fdp,DS_result$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('MultipleDataSplitting',i,DS_result$MDS_fdp,DS_result$MDS_power))
  
  knockoff_result <- knockoff(X, y, q=0.1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('Knockoff',i,knockoff_result$fdp,knockoff_result$power))
  
  BH_result <- MBHq(X, y, q=0.1, num_split)
  ResultsDataFrame=rbind(ResultsDataFrame,c('BH',i,BH_result$fdp,BH_result$power))
  
  ### save data
  return(ResultsDataFrame)}

Results=data.frame()
for(s in 1:25){
  for(i in seq(from=5,to=13,by=1)){
    Results=rbind(Results,Compare_SignalStrenght(i,s))
    print(Results)
  }
}
names(Results)

names(Results)=c('Method','SignalStrength', 'FDR','Power')

write.xlsx(Results,file='VanillaResults.xlsx')


################################use test percentage of 50############################################
Compare_SignalStrenght=function(i,s){
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
  g=ApplyTriangleLinRegTrain(X=as.data.frame(X), y, q=0.1,num_split=10, signal_index=signal_index, amountTrain=0.5, myseed = 1)
  
  ResultsDataFrame=c('LinReg DS',i, as.numeric(g$DS_fdp),as.numeric(g$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('LinReg MS',i, as.numeric(g$MDS_fdp),as.numeric(g$MDS_power)))
  
  ### Competition
  DS_result <- DS(X,y, num_split=10, q=0.1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('DataSplitting',i,DS_result$DS_fdp,DS_result$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('MultipleDataSplitting',i,DS_result$MDS_fdp,DS_result$MDS_power))
  
  knockoff_result <- knockoff(X, y, q=0.1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('Knockoff',i,knockoff_result$fdp,knockoff_result$power))
  
  BH_result <- MBHq(X, y, q=0.1, num_split)
  ResultsDataFrame=rbind(ResultsDataFrame,c('BH',i,BH_result$fdp,BH_result$power))
  
  ### save data
  return(ResultsDataFrame)}

Results=data.frame()
for(s in 1:25){
  for(i in seq(from=5,to=13,by=1)){
    Results=rbind(Results,Compare_SignalStrenght(i,s))
    print(Results)
  }
}
names(Results)

names(Results)=c('Method','SignalStrength', 'FDR','Power')

write.xlsx(Results,file='VanillaResults_train50.xlsx')


################################smaller signal strength############################################
Compare_SignalStrenght=function(i,s){
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
  g=ApplyTriangleLinRegTrain(X=as.data.frame(X), y, q=0.1,num_split=10, signal_index=signal_index, amountTrain=0.5, myseed = 1)
  
  ResultsDataFrame=c('LinReg DS',i, as.numeric(g$DS_fdp),as.numeric(g$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('LinReg MS',i, as.numeric(g$MDS_fdp),as.numeric(g$MDS_power)))
  
  ### Competition
  DS_result <- DS(X,y, num_split=10, q=0.1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('DataSplitting',i,DS_result$DS_fdp,DS_result$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('MultipleDataSplitting',i,DS_result$MDS_fdp,DS_result$MDS_power))
  
  knockoff_result <- knockoff(X, y, q=0.1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('Knockoff',i,knockoff_result$fdp,knockoff_result$power))
  
  BH_result <- MBHq(X, y, q=0.1, num_split)
  ResultsDataFrame=rbind(ResultsDataFrame,c('BH',i,BH_result$fdp,BH_result$power))
  
  ### save data
  return(ResultsDataFrame)}

Results=data.frame()
for(s in 1:25){
  for(i in seq(from=3,to=7,by=1)){
    Results=rbind(Results,Compare_SignalStrenght(i,s))
    print(Results)
  }
}
names(Results)

names(Results)=c('Method','SignalStrength', 'FDR','Power')

write.xlsx(Results,file='VanillaResults_train50.xlsx')
#Results=read.xlsx('Vanilla.xlsx')

####Visualise the results
# Create a color palette based on the number of unique methods

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
    Method == "BH" ~ "Benjamini–Hochberg",
    Method == "DataSplitting" ~ "Dai (single split)",
    Method == "MultipleDataSplitting" ~ "Dai (10 splits)",
    Method == "Knockoff" ~ "Knockoff",
    Method == "LinReg DS" ~ "Delporte (single split)",
    Method == "LinReg MS" ~ "Delporte (10 splits)",
    TRUE ~ Method  # default if none match
  ))
PowerPlot <- ggplot(resultsagg, aes(x = Signal_noisy, y = as.numeric(Avg_Power), color = Method)) +
  geom_point(size = 3) +
  geom_line()+
  labs(x = "Signal", y = "Power") +
  scale_x_continuous(breaks = seq(from = 5, to = 13, by = 1)) +
  geom_hline(yintercept = 0.8) +
  scale_color_manual(values = colors)+
  coord_cartesian(ylim = c(0.1, 1)) 
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
PlotPermute


PowerPlot <- ggplot(Results2[1:66,], aes(x = Signal_noisy, y = as.numeric(Power), color = Method)) +
  geom_point(size = 3) +
  geom_line()+
  labs(x = "Signal", y = "Power") +
  scale_x_continuous(breaks = seq(from = 3, to = 13, by = 1)) +
  geom_hline(yintercept = 0.8) +
  scale_color_manual(values = colors)
FDRPlot=ggplot(Results2[1:66,], aes(x = Signal_noisy, y = as.numeric(FDR ), color = Method)) +
  geom_point(size = 3) +
  geom_line()+
  labs(x = "Signal", y = "FDP")+
  scale_x_continuous(breaks=seq(from=3,to=13,by=1))+
  geom_hline(yintercept=0.1)+
  scale_color_manual(values = colors)
PlotPermute=ggarrange(
  PowerPlot, FDRPlot,
  common.legend = TRUE, legend = "right"
)
PlotPermute


##########apply linear boosrter###############
Compare_SignalStrenght=function(i,s){
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
  g=ApplyTriangleBoostTrain(X, y, q=0.1,num_split=10,mybooster='gblinear', signal_index=signal_index, amountTrain=0.333, myseed = 1)
  
  ResultsDataFrame=c('LinBooster DS',i, as.numeric(g$DS_fdp),as.numeric(g$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('LinBooster MS',i, as.numeric(g$MDS_fdp),as.numeric(g$MDS_power)))
  
  ### Competition
  DS_result <- DS(X,y, num_split=10, q=0.1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('DataSplitting',i,DS_result$DS_fdp,DS_result$DS_power))
  ResultsDataFrame=rbind(ResultsDataFrame,c('MultipleDataSplitting',i,DS_result$MDS_fdp,DS_result$MDS_power))
  
  knockoff_result <- knockoff(X, y, q=0.1)
  ResultsDataFrame=rbind(ResultsDataFrame,c('Knockoff',i,knockoff_result$fdp,knockoff_result$power))
  
  BH_result <- MBHq(X, y, q=0.1, num_split)
  ResultsDataFrame=rbind(ResultsDataFrame,c('BH',i,BH_result$fdp,BH_result$power))
  
  ### save data
  return(ResultsDataFrame)}

Results=data.frame()
for(s in 1:25){
  for(i in seq(from=5,to=13,by=1)){
    Results=rbind(Results,Compare_SignalStrenght(i,s))
    print(Results)
  }
}
names(Results)

names(Results)=c('Method','SignalStrength', 'FDR','Power')

#write.xlsx(Results,file='VanillaResultsBooster.xlsx')

Results=read.xlsx('VanillaResultsBooster.xlsx')
colors <- c("#000000","#FF00FF","#009900", "#99ccff", "#0000FF", "#FF0000")
Results2=Results
Results2$FDR=round(as.numeric(Results2$FDR),3)
Results2$Power=round(as.numeric(Results2$Power),2)
resultsagg <- Results2 %>%
  group_by(Method, SignalStrength) %>%
  summarize(
    Avg_FDR = mean(FDR),
    Avg_Power = mean(Power)
  )
View(resultsagg)
resultsagg$Signal_noisy <- as.numeric(resultsagg$SignalStrength) + runif(nrow(resultsagg), -0.2, 0.2)

PowerPlot <- ggplot(resultsagg, aes(x = Signal_noisy, y = as.numeric(Avg_Power), color = Method)) +
  geom_point(size = 3) +
  geom_line()+
  labs(x = "Signal", y = "Power") +
  scale_x_continuous(breaks = seq(from = 5, to = 13, by = 1)) +
  geom_hline(yintercept = 0.8) +
  scale_color_manual(values = colors)
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
PlotPermute
