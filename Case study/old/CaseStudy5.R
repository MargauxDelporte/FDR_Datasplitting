#https://waldronlab.io/curatedMetagenomicData/articles/available-studies.html
#source: https://www.bioconductor.org/packages/release/data/experiment/vignettes/curatedMetagenomicData/inst/doc/curatedMetagenomicData.html#:~:text=Abstract,collected%20from%20different%20body%20sites.
library(DT)
library(stringr)
library(mia)
library(ggplot2)
library(dplyr)
library(caret)


# Load the packages to confirm installation
library(mia)
library(curatedMetagenomicData)

# Filter metadata
Study <- sampleMetadata |>
  filter(age >= 18,
         !is.na(BMI),
         body_site == "stool",
antibiotics_current_use == "no",
gender=='female'
)|>
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
dim(final_df) #4694 1342
head(final_df[, 1:10])
names(final_df)

#keep only relevant columns
final_df$non_westernized
final_df$antibiotics_current_use
table(final_df$antibiotics_current_use) #idk if i should correct, var=4
which(grepl('species:',names(final_df)))
mydata=final_df[,c(which(names(final_df)=='age'),which(names(final_df)=='BMI'),which(grepl('species:',names(final_df))))]

#make the distribution of hscrp more normal
hist(mydata$BMI[mydata$BMI < 100],
     breaks = 30,              # adjust bin count
     main   = "Histogram of BMI (< 100)",
     xlab   = "BMI",
     col    = "skyblue",
     border = "white")

hist(log(mydata$BMI[mydata$BMI < 100]),
     breaks = 30,              # adjust bin count
     main   = "Histogram of hsCRP (< 100)",
     xlab   = "hsCRP",
     col    = "skyblue",
     border = "white")
hist(log(mydata$BMI),
     breaks = 30,              # adjust bin count
     main   = "Histogram of hsCRP (< 100)",
     xlab   = "hsCRP",
     col    = "skyblue",
     border = "white")
mydata$log_BMI=log(mydata$BMI)
sum(is.na(mydata$log_BMI))

# Calculate proportion of zeros per column (excluding outcome + IDs)
zero_prop <- colMeans(mydata == 0, na.rm = TRUE)

# Keep features with <=20% zeros OR non-numeric cols you want to preserve
keep_cols <- names(zero_prop[zero_prop <= 0.2])
mydata=mydata[,keep_cols]
#set X and y
X=mydata[,-c(which(names(mydata)=='BMI'|names(mydata)=='log_BMI'))]
ncol(X) #27
names(X)
head(X[,1:10])
y=mydata$log_BMI
n=nrow(mydata)
amountTrain=0.5
amountTest=1-amountTrain
#drop cols with more than 20% zero

train_index <- sample.int(n, size = floor(amountTrain * n), replace = FALSE)
remaining_index <- setdiff(seq_len(n), train_index)
data=cbind(y,X)
# split the remaining half evenly
size_half <- floor((amountTest/2) * n)
sample_index1 <- sample(remaining_index, size = size_half, replace = FALSE)
sample_index2 <- setdiff(remaining_index, sample_index1)
dataTrain <- data[train_index, , drop = FALSE]


mars_poly= earth(
  y ~ .,
  data    = dataTrain
)

lm <- mars_poly
lm
# --- R^2 on train (not stored) / test halves (stored) ---?earth
test1=as.data.frame(X[sample_index1,])
pred1 <- predict(lm, newdata = as.data.frame(X[sample_index1,]))
pred2 <- predict(lm, newdata = X[sample_index2,])
y1 <- y[sample_index1]; y2 <- y[sample_index2]

R2orig1 <- 1 - sum((y1 - pred1)^2) / sum((y1 - mean(y1))^2)
R2orig2 <- 1 - sum((y2 - pred2)^2) / sum((y2 - mean(y2))^2)
R2orig1;R2orig2
