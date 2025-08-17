###curated metagemomics###
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mia")
BiocManager::install("curatedMetagenomicData")
library(dplyr)
library(DT)
library(stringr)
library(mia)
#library(scater)
#library(vegan)

sampleMetadata |>
  filter(study_name == "AsnicarF_2017") |>
  select(where(~ !any(is.na(.x)))) |>
  slice(1:10) |>
  select(1:10) |>
  datatable(options = list(dom = "t"), extensions = "Responsive")


alcoholStudy <-
  filter(sampleMetadata, age >= 18) |>
  filter(!is.na(alcohol)) |>
  filter(body_site == "stool") |>
  select(where(~ !all(is.na(.x)))) |>
  returnSamples("relative_abundance", rownames = "short")

