# https://waldronlab.io/curatedMetagenomicDataAnalyses/articles/MLdatasets.html

BiocManager::install("waldronlab/curatedMetagenomicAnalyses")

library(curatedMetagenomicData)
library(curatedMetagenomicAnalyses)
library(dplyr)
library(phyloseq)
library(tidyverse)
# library(DESeq2)
library(pheatmap)


data("sampleMetadata")
availablediseases <- pull(sampleMetadata, study_condition) %>%
  table() %>%
  sort(decreasing = TRUE)
availablediseases

studies <- lapply(names(availablediseases), function(x){
  filter(sampleMetadata, study_condition %in% x) %>%
    pull(study_name) %>%
    unique()
})
names(studies) <- names(availablediseases)
studies <- studies[-grep("control", names(studies))] #get rid of controls
studies <- studies[sapply(studies, length) > 1] #available in more than one study
studies

curatedMetagenomicData("QinN_2014")

# 1. Load the dataset
gf <- curatedMetagenomicData("QinN_2014.gene_families", dryrun = FALSE)[[1]]
ma <- curatedMetagenomicData("QinN_2014.marker_abundance", dryrun = FALSE)[[1]]
mp <- curatedMetagenomicData("QinN_2014.marker_presence", dryrun = FALSE)[[1]]
pa <- curatedMetagenomicData("QinN_2014.pathway_abundance", dryrun = FALSE)[[1]]
pc <- curatedMetagenomicData("QinN_2014.pathway_coverage", dryrun = FALSE)[[1]]
ra_list <- curatedMetagenomicData("QinN_2014.relative_abundance", dryrun = FALSE)
ra_list
ra <- ra_list[[1]]   # your line is fine
class(ra)
# 2. Look at the metadata (phenotype data)
pd <- as.data.frame(colData(ra))
head(pd)
str(pd)
pd$sample_id <- rownames(pd)


#3. abundance matrix
assayNames(ra)      # see what the assay is called, e.g. "abundance", "counts", etc.

ra_mat <- assay(ra) # default assay
dim(ra_mat)

#4. feature info
taxa_info <- as.data.frame(rowData(ra))
head(taxa_info)
all(colnames(ra_mat) == rownames(pd))
