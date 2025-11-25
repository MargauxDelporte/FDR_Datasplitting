## ----------------------------------------------------------
## Packages
## ----------------------------------------------------------
library(curatedMetagenomicData)
library(TreeSummarizedExperiment)
library(tidyverse)
library(vegan)
library(pheatmap)

## ----------------------------------------------------------
## Load QinN_2014 relative abundance (TreeSummarizedExperiment)
## ----------------------------------------------------------
ra_list <- curatedMetagenomicData("QinN_2014.relative_abundance", dryrun = FALSE)
ra <- ra_list[[1]]   # 2021-03-31.QinN_2014.relative_abundance

## ----------------------------------------------------------
## Extract sample metadata, abundance matrix, taxonomy
## ----------------------------------------------------------
pd <- as.data.frame(colData(ra))
pd$sample_id <- rownames(pd)

ra_mat <- assay(ra)      # taxa x samples
tax <- as.data.frame(rowData(ra))
tax$feature_id <- rownames(tax)

## Inspect available metadata variables
str(pd)
names(pd)

## ----------------------------------------------------------
## Choose outcome variable (assume 'disease' exists)
## ----------------------------------------------------------
# If name is different, change "disease" here:
outcome_var <- "disease"

pd <- pd %>%
  filter(!is.na(.data[[outcome_var]]))

# Keep only those samples in abundance matrix
ra_mat <- ra_mat[, pd$sample_id]

pd[[outcome_var]] <- factor(pd[[outcome_var]])
table(pd[[outcome_var]])

## ----------------------------------------------------------
## Alpha diversity (Shannon + observed richness)
## ----------------------------------------------------------
shannon <- diversity(t(ra_mat), index = "shannon")
observed <- colSums(ra_mat > 0)

alpha_df <- data.frame(
  sample_id = colnames(ra_mat),
  shannon   = shannon,
  observed  = observed
) %>%
  left_join(pd, by = "sample_id")

# Shannon
ggplot(alpha_df, aes(x = .data[[outcome_var]], y = shannon,
                     fill = .data[[outcome_var]])) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) +
  theme_bw() +
  labs(x = outcome_var, y = "Shannon diversity",
       title = "Shannon diversity – QinN_2014")

# Observed richness
ggplot(alpha_df, aes(x = .data[[outcome_var]], y = observed,
                     fill = .data[[outcome_var]])) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) +
  theme_bw() +
  labs(x = outcome_var, y = "Observed taxa",
       title = "Observed richness – QinN_2014")

## ----------------------------------------------------------
## Beta diversity (Bray–Curtis) + MDS plot + PERMANOVA
## ----------------------------------------------------------
dist_bc <- vegdist(t(ra_mat), method = "bray")

ord <- cmdscale(dist_bc, k = 2)
ord_df <- data.frame(
  sample_id = rownames(ord),
  MDS1 = ord[, 1],
  MDS2 = ord[, 2]
) %>%
  left_join(pd, by = "sample_id")

ggplot(ord_df, aes(x = MDS1, y = MDS2,
                   color = .data[[outcome_var]])) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "Bray–Curtis MDS – QinN_2014")

# PERMANOVA
adonis_res <- adonis2(dist_bc ~ pd[[outcome_var]])
adonis_res

## ----------------------------------------------------------
## Collapse to species level and explore composition
## ----------------------------------------------------------
# Row names look like "k__Bacteria|p__...|s__Species_name"
tax_strings <- rownames(ra_mat)

species <- ifelse(
  grepl("s__", tax_strings),
  sub(".*\\|s__", "", tax_strings),
  "Unclassified"
)

# Collapse abundances by species
species_abun <- rowsum(ra_mat, group = species)

# Relative abundances at species level (normalize columns)
species_rel <- sweep(species_abun, 2, colSums(species_abun), "/")

# Top 20 most abundant species
top_species <- names(sort(rowMeans(species_rel), decreasing = TRUE))[1:20]

species_long <- species_rel[top_species, ] %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  pivot_longer(-species, names_to = "sample_id", values_to = "rel_abun") %>%
  left_join(pd, by = "sample_id")

ggplot(species_long,
       aes(x = sample_id, y = rel_abun, fill = species)) +
  geom_col() +
  facet_grid(rows = vars(.data[[outcome_var]]), scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Samples", y = "Relative abundance",
       title = "Top 20 species – QinN_2014")

## ----------------------------------------------------------
## Heatmap of top variable species
## ----------------------------------------------------------
# Variance by species
species_var <- apply(species_rel, 1, var)
hv_species <- names(sort(species_var, decreasing = TRUE))[1:30]

anno_col <- data.frame(group = pd[[outcome_var]])
rownames(anno_col) <- pd$sample_id

pheatmap(
  species_rel[hv_species, pd$sample_id],
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",  # <- changed from "bray"
  annotation_col = anno_col,
  main = "Top variable species – QinN_2014"
)

## ----------------------------------------------------------
## Simple differential abundance (Wilcoxon test) at species level
## (Assumes 2-level outcome)
## ----------------------------------------------------------?anova
group <- pd[[outcome_var]]
group <- droplevels(group)
stopifnot(nlevels(group) == 2)
lvl <- levels(group)

da_list <- lapply(rownames(species_abun), function(sp) {
  x <- species_abun[sp, pd$sample_id]
  # add small pseudocount for log2FC
  w <- lm(x ~ group)
  data.frame(
    species = sp,
    p       = w$p.value,
    log2FC  = log2(mean(x[group == lvl[2]] + 1e-6) /
                     mean(x[group == lvl[1]] + 1e-6))
  )
})

da_res <- bind_rows(da_list) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH")
  ) %>%
  arrange(p_adj)

head(da_res, 20)

# Volcano plot
ggplot(da_res, aes(x = log2FC, y = -log10(p_adj))) +
  geom_point(alpha = 0.4) +
  theme_bw() +
  labs(
    x = paste0("log2( ", lvl[2], " / ", lvl[1], " )"),
    y = "-log10(FDR)",
    title = "Differential abundance (Wilcoxon, species level)"
  )
