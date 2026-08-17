##volgende test: wat ophogen sample size en zien of het dan hogere f aankan

###Linear model
rm(list = ls())
#install.packages("remotes")  # if needed
#remotes::install_github("cran/hdi")
#mywd='/home/mde4023/FDR_Datasplitting'
#mywd='C:/Users/mde4023/Downloads/FDR_Datasplitting'
mywd='C:/Users/mraga/Downloads/FDR_Datasplitting'
setwd(mywd)
source(paste0(mywd,'/Functions/HelperFunctions.R'))
source(paste0(mywd,'/Scenario/Scenario5/TriangleRidgeRegTrainMS.R'))
source(paste0(mywd,'/Functions Dai/knockoff.R'))
source(paste0(mywd,'/Functions Dai/analysis.R'))
source(paste0(mywd,'/Functions Dai/MBHq2.R'))
source(paste0(mywd,'/Functions Dai/DS.R'))
source(paste0(mywd,'/Functions Dai/fdp_power.R'))
library(mvtnorm)
library(MASS)
library(glmnet)
library(knockoff)
library(mvtnorm)
library(parallel)
#install.packages("hdi",build = FALSE) 
### algorithmic settings
num_split <- 50
n <- 1500
p <- 250
p0 <- 25
q <- 0.1
#rho <- 0.25
set.seed(557)
signal_index <- sample(c(1:p), size = p0, replace = F)
#View(cor(X))
#s=runif(1)*10000000
#i=8

Compare_SignalStrength <- function(i, s, f) {
  set.seed(s)
  delta <- i
  # simulate data
  n1 <- floor(n/2); n2 <- n - n1
  
  # Target equicorrelation matrix
  n_df <- f*p  # much larger than p+1
  D <- matrix(rWishart(n = 1, df = n_df, Sigma = diag(p)), nrow = p)
  R <- cov2cor(D)
  
  X1 <- rmvnorm(n1, rep(0, p), sigma = D)
  X2 <- rmvnorm(n2, rep(0, p), sigma = D)
  X  <- rbind(X1, X2)
  beta_star <- rep(0, p)
  beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta * sqrt(log(p) / n))
  
  ### generate y
  y <- X %*% beta_star + rnorm(n, mean = 0, sd = 1)
  
  ### Compute VIF from empirical correlation matrix of X
  ## VIF_j = [cor(X)^{-1}]_jj
  R_emp <- cor(X)
  vif_vals <- tryCatch(
    diag(solve(R_emp)),
    error = function(e) rep(NA, p)  ## in case R_emp is singular/near-singular
  )
  mean_VIF <- round(mean(vif_vals, na.rm = TRUE),2)
  max_VIF  <- round(max(vif_vals, na.rm = TRUE),2)
  
  ###my own methods:
  g <- ApplyTriangleRidgeTrain(X = as.data.frame(X), y, q = q, num_split = num_split,
                               signal_index = signal_index, myseed = 5)
  
  ResultsDataFrame <- c('LinReg DS', i, round(as.numeric(g$DS_fdp),2),  round(as.numeric(g$DS_power),2), mean_VIF, max_VIF)
  ResultsDataFrame <- rbind(ResultsDataFrame, c('LinReg MS', i,  round(as.numeric(g$MDS_fdp),2),  round(as.numeric(g$MDS_power),2), mean_VIF, max_VIF))
  
  ### Competition
  DS_result <- DS(X, y, num_split = num_split, q = q)
  ResultsDataFrame <- rbind(ResultsDataFrame, c('DataSplitting', i, DS_result$DS_fdp, DS_result$DS_power, mean_VIF, max_VIF))
  ResultsDataFrame <- rbind(ResultsDataFrame, c('MultipleDataSplitting', i, DS_result$MDS_fdp, DS_result$MDS_power, mean_VIF, max_VIF))
  
  knockoff_result <- knockoff(X, y, q = q)
  ResultsDataFrame <- rbind(ResultsDataFrame, c('Knockoff', i, knockoff_result$fdp, knockoff_result$power, mean_VIF, max_VIF))
  
  BH_result <- MBHq(X, y, q = q, num_split = num_split)
  ResultsDataFrame <- rbind(ResultsDataFrame, c('BH', i, BH_result$fdp, BH_result$power, mean_VIF, max_VIF))
  
  ### save data
  colnames(ResultsDataFrame) <- c('Method', 'delta', 'fdp', 'power', 'mean_VIF', 'max_VIF')
  return(ResultsDataFrame)
}

Compare_SignalStrength(i=7,s=559,f=5)

# Load required packages
pkgs <- c('MASS','glmnet','knockoff','mvtnorm',
          'foreach','doParallel')
lapply(pkgs, library, character.only = TRUE)

# === PARAMETER GRID ===
param_grid <- expand.grid(
  s = 17:200,
  i = 13
)

ncore=30
# === SET UP PARALLEL BACKEND ===
cl <- makeCluster(ncore-1)
# export working dir so workers can source
clusterExport(cl, 'mywd')
# have each worker source & load libraries
clusterEvalQ(cl, {
  setwd(mywd)
  source(file.path(mywd, '/Functions', 'HelperFunctions.R'))
  source(paste0(mywd,'/Scenario/Scenario5/TriangleRidgeRegTrainMS.R'))
  source(file.path(mywd, '/Functions Dai', 'knockoff.R'))
  source(file.path(mywd, '/Functions Dai', 'analysis.R'))
  source(file.path(mywd, '/Functions Dai', 'MBHq2.R'))
  source(file.path(mywd, '/Functions Dai', 'DS.R'))
  source(file.path(mywd, '/Functions Dai', 'fdp_power.R'))
  lapply( c('MASS','glmnet','knockoff','mvtnorm',
            'foreach','doParallel'),
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
  chunk <- Compare_SignalStrength(i = i_val, s = s_val,f=5)
  
  # write out this chunk immediately
  fname <- sprintf("Results_s%02d_i%02d.csv", s_val, i_val)
  write.csv(chunk, file = paste0(mywd,"/Scenario/Scenario5/Temp2/",fname), row.names = FALSE)
  
  # return for final binding
  chunk
}
# === CLEANUP AND FINAL SAVE ===
stopCluster(cl)
warnings()


library(openxlsx)
library(dplyr)

# Set your working directory
mywd2 <- paste0(mywd,'/Scenario/Scenario5/Temp2/')

# Get list of all Excel files (.xlsx and .xls)
excel_files <- list.files(path = mywd2, pattern = "\\.csv?$", full.names = TRUE)

# Read and combine all files
combined_data <- excel_files %>%
  lapply(read.csv) %>%
  bind_rows()
names(combined_data)[2:4]=c('SignalStrength','FDR','Power')
combined_data$Method <- dplyr::recode(
  combined_data$Method,
  "LinReg DS"             = "Delporte (single split)",
  "LinReg MS"             = "Delporte (50 splits)",
  "DataSplitting"         = "Dai (single split)",
  "MultipleDataSplitting" = "Dai (50 splits)",
  "Knockoff"              = "Knockoff",
  "BH"                    = "Benjamini–Hochberg (BH)"
)

write.xlsx(combined_data,paste0(mywd,'/Scenario/Scenario5/','Scenario5.xlsx'))
s5=combined_data


head(s5)
resultsagg <- s5 %>%
  group_by(Method, SignalStrength) %>%
  summarize(
    Avg_FDR = mean(FDR, na.rm = TRUE),
    Avg_Power = mean(Power, na.rm = TRUE),
    # Add empirical confidence intervals
    FDR_SE = sd(FDR, na.rm = TRUE) / sqrt(sum(!is.na(FDR))),
    Power_SE = sd(Power, na.rm = TRUE) / sqrt(sum(!is.na(Power))),
    FDR_Lower = Avg_FDR - 1.96 * FDR_SE,
    FDR_Upper = Avg_FDR + 1.96 * FDR_SE,
    Power_Lower = Avg_Power - 1.96 * Power_SE,
    Power_Upper = Avg_Power + 1.96 * Power_SE,
    # Alternative: use quantiles for empirical CI
    FDR_Q025 = quantile(FDR, 0.025, na.rm = TRUE),
    FDR_Q975 = quantile(FDR, 0.975, na.rm = TRUE),
    Power_Q025 = quantile(Power, 0.025, na.rm = TRUE),
    Power_Q975 = quantile(Power, 0.975, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  ) %>%
  mutate(Scenario = 'scenario_name')  # Add scenario label

print(resultsagg)
View(resultsagg)
