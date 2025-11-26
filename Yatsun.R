# 1. Install the package if you haven't
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("curatedMetagenomicData")
BiocManager::install("phyloseq")

# 2. Load the library and fetch the dataset
library(curatedMetagenomicData)
library(phyloseq)

# Fetch the specific dataset (YatsunenkoT_2012)
# dryrun=FALSE downloads the actual data
dataset_list <- curatedMetagenomicData("YatsunenkoT_2012.relative_abundance", dryrun = FALSE)
curatedMetagenomicData("Yatsunenko", dryrun = TRUE)
# 3. Extract the actual data object
tse <- dataset_list[[1]]

# 4. Access the features (Taxa) and Outcome (Age)
# The abundance table (Features)
abundance_table <- assays(tse)[["relative_abundance"]]

# The metadata (Outcome: age)
metadata <- colData(tse)
continuous_outcome <- metadata$age

# Check dimensions
print(dim(abundance_table)) 
# Returns approx: [Rows=Taxa, Cols=Samples]