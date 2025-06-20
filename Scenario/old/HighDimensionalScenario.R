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
source(paste0(mywd,'/Functions Dai/FDR_power.R'))

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
    FDR    = numeric(),
    Power  = numeric(),
    stringsAsFactors = FALSE
  )
  
  # bind all rows
  ResultsDataFrame <- rbind(
    ResultsDataFrame,
    data.frame(seed=s, Method = "Boost DS",                Delta = i, FDR = g1$DS_FDR,    Power = g1$DS_power),
    data.frame(seed=s, Method = "Boost MS",                Delta = i, FDR = g1$MDS_FDR,   Power = g1$MDS_power),
    data.frame(seed=s, Method = "DataSplitting",           Delta = i, FDR = DS_result$DS_FDR,  Power = DS_result$DS_power),
    data.frame(seed=s, Method = "MultipleDataSplitting",   Delta = i, FDR = DS_result$MDS_FDR, Power = DS_result$MDS_power),
    data.frame(seed=s, Method = "Knockoff",                Delta = i, FDR = knockoff_result$FDR, Power = knockoff_result$power),
    data.frame(seed=s, Method = "Benjamini–Hochberg (BH)", Delta = i, FDR = BH_result$FDR,     Power = BH_result$power)
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
source(file.path(mywd, 'Functions Dai', 'FDR_power.R'))

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
  source(file.path(mywd, 'Functions Dai', 'FDR_power.R'))
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

##########read in the results#############
data_dir <- file.path(mywd, "Temp")

# initialize a list to collect each chunk
results_list <- vector("list", length = (13-7+1) * 50)
idx <- 1

for (s in 1:50) {
  for (i in 7:13) {
    # build filename and full path
    fname <- sprintf("Results_s%02d_i%02d.csv", s, i)
    fpath <- file.path(data_dir, fname)
    
    # error out if the file doesn't exist
    if (!file.exists(fpath)) {
      stop(sprintf("ERROR: File not found: %s", fpath))
    }
    
    # read it in
    df <- read.csv(fpath, stringsAsFactors = FALSE)
    
    # tag with s and i
    df$s <- s
    df$i <- i
    
    # store into our list
    results_list[[idx]] <- df
    idx <- idx + 1
  }
}

# combine all into one data.frame
all_results <- do.call(rbind, results_list)

# inspect
head(all_results)

##########visualise the results###########
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readr)
library(openxlsx)
write.xlsx(all_results,'High_Dimensional_results.xlsx')
#Results=read.xlsx('VanillaResults.xlsx')
colors <- c("#000000","#FF00FF","#009900", "#99ccff", "#0000FF", "#FF0000")
Results2=all_results
names(Results2)=c('seed','Method','SignalStrength','FDR','Power','s','i')
Results2$FDR=round(as.numeric(Results2$FDR),3)
Results2$Power=round(as.numeric(Results2$Power),2)
resultsagg <- Results2 %>%
  group_by(Method, SignalStrength) %>%
  summarize(
    Avg_FDR = mean(FDR),
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
  labs(x = "Signal", y = "FDR")+
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
