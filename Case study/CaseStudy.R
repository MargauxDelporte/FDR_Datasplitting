#source: https://www.bioconductor.org/packages/release/data/experiment/vignettes/curatedMetagenomicData/inst/doc/curatedMetagenomicData.html#:~:text=Abstract,collected%20from%20different%20body%20sites.
library(dplyr)
library(DT)
library(stringr)
library(mia)
library(scater)
library(vegan)
library(ggplot2)

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mia", type = "source")
# Update all old packages without asking
update.packages(ask = FALSE)

# Install 'mia' from Bioconductor
if (!requireNamespace("mia", quietly = TRUE))
  BiocManager::install("mia", ask = FALSE, update = TRUE)

# Install 'curatedMetagenomicData' from Bioconductor
if (!requireNamespace("curatedMetagenomicData", quietly = TRUE))
  BiocManager::install("curatedMetagenomicData", ask = FALSE, update = TRUE)

# Load the packages to confirm installation
library(mia)
library(curatedMetagenomicData)
names(sampleMetadata)
tab <- table(sampleMetadata$disease)
tab_filtered <- tab[tab > 200] 
tab_filtered
table(sampleMetadata$disease)
Study <-
  filter(sampleMetadata, age >= 18) |>
  filter(!is.na(hscrp)) |>
  filter(body_site == "stool") |>
  select(where(~ !all(is.na(.x)))) |>
  returnSamples("relative_abundance", rownames = "short")

altExps(Study) <-
  splitByRanks(Study)

Study |>
  estimateDiversity(assay.type = "relative_abundance", index = "shannon") |>
  plotColData(x = "alcohol", y = "shannon", colour_by = "alcohol", shape_by = "alcohol") +
  labs(x = "Alcohol", y = "Alpha Diversity (H')") +
  guides(colour = guide_legend(title = "Alcohol"), shape = guide_legend(title = "Alcohol")) +
  theme(legend.position = "none")
Study |>
  runMDS(FUN = vegdist, method = "bray", exprs_values = "relative_abundance", altexp = "genus", name = "BrayCurtis") |>
  plotReducedDim("BrayCurtis", colour_by = "alcohol", shape_by = "alcohol") +
  labs(x = "PCo 1", y = "PCo 2") +
  guides(colour = guide_legend(title = "Alcohol"), shape = guide_legend(title = "Alcohol")) +
  theme(legend.position = c(0.90, 0.85))
Study |>
  runUMAP(exprs_values = "relative_abundance", altexp = "genus", name = "UMAP") |>
  plotReducedDim("UMAP", colour_by = "alcohol", shape_by = "alcohol") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(colour = guide_legend(title = "Alcohol"), shape = guide_legend(title = "Alcohol")) +
  theme(legend.position = c(0.90, 0.85))
