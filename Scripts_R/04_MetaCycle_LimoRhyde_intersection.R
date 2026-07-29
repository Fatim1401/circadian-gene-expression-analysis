
############################################################
# Script: 04_MetaCycle_LimoRhyde_intersection.R
# Project: Circadian transcriptomic analysis under abiotic stress
#
# Description:
# Comparison of rhythmic genes identified by MetaCycle and
# LimoRhyde approaches. The analysis identifies overlapping
# rhythmic genes detected by both methods in leaves and seeds.
#
# Input:
# - MetaCycle rhythmic gene lists
# - LimoRhyde rhythmic gene lists
#
# Output:
# - Common rhythmic genes between methods
# - Intersection tables for seeds and leaves
#
# Author: Fatoumata KAMISSOKO 
############################################################



## Load libraries -------------------------------------------------------------

library(data.table)



## Define paths ---------------------------------------------------------------

data_dir <- "outputs/data"



## Load MetaCycle results -----------------------------------------------------


meta_seed <- read.table(
  file.path(
    data_dir,
    "Cycling_seeds.txt"
  ),
  header = TRUE,
  sep = "\t"
)


meta_leaf <- read.table(
  file.path(
    data_dir,
    "Cycling_leaf.txt"
  ),
  header = TRUE,
  sep = "\t"
)



## Load LimoRhyde results -----------------------------------------------------


limo_seed <- read.table(
  file.path(
    data_dir,
    "LimoRhyde_rhythmic_genes_seeds.txt"
  ),
  header = TRUE,
  sep = "\t"
)



limo_leaf <- read.table(
  file.path(
    data_dir,
    "LimoRhyde_rhythmic_genes_leaves.txt"
  ),
  header = TRUE,
  sep = "\t"
)



############################################################
# 1. Intersection in seeds
############################################################


# Extract gene identifiers

meta_seed_genes <- meta_seed$CycID


limo_seed_genes <- limo_seed$gene_id



# Common rhythmic genes

overlap_seed <- intersect(
  meta_seed_genes,
  limo_seed_genes
)



cat(
  "Common rhythmic genes in seeds:",
  length(overlap_seed),
  "\n"
)



# Create table

overlap_seed_df <- data.frame(
  gene_id = overlap_seed,
  Organ = "Seeds"
)



# Add MetaCycle phase information

overlap_seed_df <- merge(
  overlap_seed_df,
  meta_seed[, c(
    "CycID",
    "JTK_adjphase",
    "JTK_period"
  )],
  by.x = "gene_id",
  by.y = "CycID",
  all.x = TRUE
)



############################################################
# 2. Intersection in leaves
############################################################


# Extract gene identifiers

meta_leaf_genes <- meta_leaf$CycID


limo_leaf_genes <- limo_leaf$gene_id



# Common rhythmic genes

overlap_leaf <- intersect(
  meta_leaf_genes,
  limo_leaf_genes
)



cat(
  "Common rhythmic genes in leaves:",
  length(overlap_leaf),
  "\n"
)



# Create table

overlap_leaf_df <- data.frame(
  gene_id = overlap_leaf,
  Organ = "Leaves"
)



# Add MetaCycle phase information

overlap_leaf_df <- merge(
  overlap_leaf_df,
  meta_leaf[, c(
    "CycID",
    "JTK_adjphase",
    "JTK_period"
  )],
  by.x = "gene_id",
  by.y = "CycID",
  all.x = TRUE
)



############################################################
# Export results
############################################################


write.table(
  overlap_seed_df,
  file.path(
    data_dir,
    "MetaCycle_LimoRhyde_overlap_seeds.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



write.table(
  overlap_leaf_df,
  file.path(
    data_dir,
    "MetaCycle_LimoRhyde_overlap_leaves.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



# Export gene lists only

write.table(
  overlap_seed,
  file.path(
    data_dir,
    "common_rhythmic_genes_seeds.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)



write.table(
  overlap_leaf,
  file.path(
    data_dir,
    "common_rhythmic_genes_leaves.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)



############################################################
# Summary
############################################################


cat(
  "---------------------------------\n"
)

cat(
  "Seeds:\n"
)

cat(
  "MetaCycle rhythmic genes:",
  length(meta_seed_genes),
  "\n"
)

cat(
  "LimoRhyde rhythmic genes:",
  length(limo_seed_genes),
  "\n"
)

cat(
  "Common genes:",
  length(overlap_seed),
  "\n\n"
)



cat(
  "Leaves:\n"
)

cat(
  "MetaCycle rhythmic genes:",
  length(meta_leaf_genes),
  "\n"
)

cat(
  "LimoRhyde rhythmic genes:",
  length(limo_leaf_genes),
  "\n"
)

cat(
  "Common genes:",
  length(overlap_leaf),
  "\n"
)