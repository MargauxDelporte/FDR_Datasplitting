#source: https://www.bioconductor.org/packages/release/data/experiment/vignettes/curatedMetagenomicData/inst/doc/curatedMetagenomicData.html#:~:text=Abstract,collected%20from%20different%20body%20sites.
library(dplyr)
library(DT)
library(stringr)
library(mia)
library(scater)
library(vegan)
library(ggplot2)
library(SummarizedExperiment)
library(dplyr)
library(caret)
# Install BiocManager if not present
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("mia", type = "source")
# Update all old packages without asking
#update.packages(ask = FALSE)

# Install 'mia' from Bioconductor
#if (!requireNamespace("mia", quietly = TRUE))
#  BiocManager::install("mia", ask = FALSE, update = TRUE)

# Install 'curatedMetagenomicData' from Bioconductor
#if (!requireNamespace("curatedMetagenomicData", quietly = TRUE))
#  BiocManager::install("curatedMetagenomicData", ask = FALSE, update = TRUE)

# Load the packages to confirm installation
library(mia)
library(curatedMetagenomicData)

# Filter metadata
Study <- sampleMetadata |>
  filter(age >= 18,
         !is.na(hscrp),
         body_site == "stool") |>
  select(where(~ !all(is.na(.x))))

# Attach microbiome abundances
Study_with_microbiome <- returnSamples(
  Study,
  dataType = "relative_abundance",
  rownames = "short"
)

# Extract microbiome abundances (taxa x samples matrix)
microbiome_df <- as.data.frame(t(assay(Study_with_microbiome)))
microbiome_df$sample_id <- rownames(microbiome_df)

# Extract metadata (including hscrp)
metadata_df <- as.data.frame(colData(Study_with_microbiome))
metadata_df$sample_id <- rownames(metadata_df)

# Merge into one dataframe
final_df <- dplyr::left_join(metadata_df, microbiome_df, by = "sample_id")
dim(final_df)
head(final_df[, 1:10])
names(final_df)
#keep only relevant columns
final_df$non_westernized
final_df$antibiotics_current_use
table(final_df$antibiotics_current_use) #idk if i should correct, var=4

names(final_df)
mydata=final_df[,c(which(names(final_df)=='hscrp'),which(names(final_df)=='age'),which(names(final_df)=='BMI'),64:847)]
head(mydata[,1:10])
View(final_df)
#make the distribution of hscrp more normal
hist(mydata$hscrp[mydata$hscrp < 100],
     breaks = 30,              # adjust bin count
     main   = "Histogram of hsCRP (< 100)",
     xlab   = "hsCRP",
     col    = "skyblue",
     border = "white")

hist(log(mydata$hscrp[mydata$hscrp < 100]),
     breaks = 30,              # adjust bin count
     main   = "Histogram of hsCRP (< 100)",
     xlab   = "hsCRP",
     col    = "skyblue",
     border = "white")
hist(log(mydata$hscrp),
     breaks = 30,              # adjust bin count
     main   = "Histogram of hsCRP (< 100)",
     xlab   = "hsCRP",
     col    = "skyblue",
     border = "white")
mydata$log_hscrp=log(mydata$hscrp)
sum(is.na(mydata$log_hscrp))

# Calculate proportion of zeros per column (excluding outcome + IDs)
zero_prop <- colMeans(mydata == 0, na.rm = TRUE)
View(mydata)
# Keep features with <=20% zeros OR non-numeric cols you want to preserve
keep_cols <- names(zero_prop[zero_prop <= 0.2])
mydata=mydata[,keep_cols]
#set X and y
X=mydata[,-1]
ncol(X)
head(X[,1:10])
y=mydata[,1]
y

#drop cols with more than 20% zero

#crossvalidate my xgboost tree model
set.seed(123)
trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[trainIndex, ]
X_test  <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test  <- y[-trainIndex]


ctrl <- trainControl(
  method = "cv",    # k-fold cross-validation
  number = 5,
  verboseIter = FALSE
)

# Define parameter grid
grid <- expand.grid(
  nrounds = 500,
  max_depth = c(3, 5, 7),
  eta = c(0.01, 0.1, 0.3),
  gamma = c(0,0.5, 1),
  colsample_bytree = 1,
  min_child_weight = c(1, 3, 5),
  subsample = c(0.8, 1)
)

set.seed(123)
xgb_model <- train(
  x = X_train,
  y = y_train,
  method = "xgbTree",
  trControl = ctrl,
  tuneGrid = grid,
  metric = "RMSE"
)

# Best parameters
xgb_model$bestTune

# Predict on test set
pred <- predict(xgb_model, newdata = X_test)

# RMSE and R²
RMSE(pred, y_test)
R2(pred, y_test)


