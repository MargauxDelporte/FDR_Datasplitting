library(readr)
library(doParallel)
library(doRNG)      # Reproducible parallel random numbers
library(microbiome)
library(earth)
library(dplyr)
library(earth)        # MARS
library(foreach)
library(doParallel)
library(doRNG)  

source(paste0('C:/Users/mde4023/Downloads/FDR_Datasplitting','/Functions/HelperFunctions.R'))
# ==============================================================================
#Helper functions
# ==============================================================================

# Function to calculate Permuted R2 for a specific feature j
permR2 <- function(data, Y, j, model) {
  Xperm <- data
  # Permute column j
  Xperm[, j] <- sample(data[, j], replace = FALSE)
  
  # Predict using permuted data
  pred_perm <- predict(model, newdata = as.data.frame(Xperm))
  
  # Calculate R2
  rsq_perm <- 1 - sum((Y - pred_perm)^2) / sum((Y - mean(Y))^2)
  return(rsq_perm)
}


# ==============================================================================
# 1. File Paths & Data Import
# ==============================================================================

base_dir <- "C:/Users/mde4023/Downloads/FDR_Datasplitting/Case study/yatsunenko"

task <- read.delim(
  file         = file.path(base_dir, "task-baby-age.txt"),
  header       = FALSE,
  comment.char = "#"
)
task2 <- read.delim(
  file         = file.path(base_dir, "task-malawi-venezuela.txt"),
  header       = FALSE,
  comment.char = "#"
)
task3 <- read.delim(
  file         = file.path(base_dir, "task-sex.txt"),
  header       = FALSE,
  comment.char = "#"
)

merge(task,task3,by='V1')
f <- read_excel("C:/Users/mde4023/Downloads/41586_2012_BFnature11053_MOESM27_ESM/nature11053-s2/2011-02-02467D-MS2011-02-02467CYatsunenko_TableS2.xls")


USApeople=f[which(f$Country=='USA'),]

names(USApeople)[6]='Age'
USApeople$Age2=as.numeric(USApeople$Age)
USApeople[which(is.na(USApeople$Age2)),]

USAInf=USApeople[which(USApeople$Age2<4),]
View(USAInf)
!USAInf[,5]%in%task$V1
task$V1
as.vector(USAInf[,5])

task_ids <- sub("\\.[0-9]+$", "", task$V1)           # remove trailing .418xxx
usa_ids  <- as.vector(USAInf[["Sample Identifier"]])

setdiff(usa_ids, task_ids)   # in USAInf but not in task
setdiff(task_ids, usa_ids)   # in task but not in USAInf
21022