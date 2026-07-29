############################################################
# Script: 06_Rhythmicity_visualization.R
# Project: Circadian transcriptomic gene expression analysis 
#
# Description:
# Visualization and summary of rhythmic gene identification.
# This script generates Venn diagrams comparing rhythmic genes
# identified by MetaCycle and LimoRhyde, and highlights genes
# showing differential rhythmicity between leaves and seeds.
#
# Input:
# - MetaCycle rhythmic gene lists
# - LimoRhyde rhythmic gene lists
# - Differential rhythmicity results
#
# Output:
# - Venn diagrams
# - Overlap gene lists
# - Summary tables
#
# Author: Fatoumata KAMISSOKO 
############################################################



## Load libraries -------------------------------------------------------------

library(VennDiagram)
library(grid)
library(data.table)



## Define paths ---------------------------------------------------------------

data_dir <- "outputs/data"
figure_dir <- "outputs/figures"



############################################################
# 1. MetaCycle vs LimoRhyde overlap
############################################################


# Load results

meta_seed <- read.table(
  file.path(data_dir, "Cycling_seeds.txt"),
  header = TRUE,
  sep = "\t"
)

meta_leaf <- read.table(
  file.path(data_dir, "Cycling_leaf.txt"),
  header = TRUE,
  sep = "\t"
)


limo_seed <- read.table(
  file.path(data_dir,
            "LimoRhyde_rhythmic_genes_seeds.txt"),
  header = TRUE,
  sep = "\t"
)


limo_leaf <- read.table(
  file.path(data_dir,
            "LimoRhyde_rhythmic_genes_leaves.txt"),
  header = TRUE,
  sep = "\t"
)



# Extract gene lists

meta_seed_genes <- meta_seed$CycID
meta_leaf_genes <- meta_leaf$CycID

limo_seed_genes <- limo_seed$gene_id
limo_leaf_genes <- limo_leaf$gene_id



# Calculate overlap

overlap_seed <- intersect(
  meta_seed_genes,
  limo_seed_genes
)


overlap_leaf <- intersect(
  meta_leaf_genes,
  limo_leaf_genes
)



############################################################
# Venn diagram seeds
############################################################


cairo_pdf(
  file.path(
    figure_dir,
    "Venn_MetaCycle_LimoRhyde_seeds.pdf"
  ),
  width = 6,
  height = 6
)


draw.pairwise.venn(
  area1 = length(meta_seed_genes),
  area2 = length(limo_seed_genes),
  cross.area = length(overlap_seed),
  category = c(
    "MetaCycle",
    "LimoRhyde"
  ),
  alpha = 0.7,
  lty = "blank"
)


dev.off()



############################################################
# Venn diagram leaves
############################################################


cairo_pdf(
  file.path(
    figure_dir,
    "Venn_MetaCycle_LimoRhyde_leaves.pdf"
  ),
  width = 6,
  height = 6
)


draw.pairwise.venn(
  area1 = length(meta_leaf_genes),
  area2 = length(limo_leaf_genes),
  cross.area = length(overlap_leaf),
  category = c(
    "MetaCycle",
    "LimoRhyde"
  ),
  alpha = 0.7,
  lty = "blank"
)


dev.off()



############################################################
# 2. Differential rhythmicity summary
############################################################


dr <- read.table(
  file.path(
    data_dir,
    "differentially_rhythmic_genes.txt"
  ),
  header = TRUE,
  sep = "\t"
)



# Common rhythmic genes between organs

common_rhythmic <- intersect(
  overlap_seed,
  overlap_leaf
)



common_DR <- intersect(
  common_rhythmic,
  dr$gene_id
)



cat(
  "Common rhythmic genes:",
  length(common_rhythmic),
  "\n"
)


cat(
  "Common differentially rhythmic genes:",
  length(common_DR),
  "\n"
)



############################################################
# Export final gene lists
############################################################


write.table(
  data.frame(
    gene_id = common_rhythmic
  ),
  file.path(
    data_dir,
    "common_rhythmic_leaf_seed.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



write.table(
  data.frame(
    gene_id = common_DR
  ),
  file.path(
    data_dir,
    "common_differentially_rhythmic_genes.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



############################################################
# Final summary
############################################################


summary_table <- data.frame(
  Category = c(
    "MetaCycle seeds",
    "MetaCycle leaves",
    "LimoRhyde seeds",
    "LimoRhyde leaves",
    "Common seeds",
    "Common leaves",
    "Differential rhythmic genes"
  ),
  
  Number_of_genes = c(
    length(meta_seed_genes),
    length(meta_leaf_genes),
    length(limo_seed_genes),
    length(limo_leaf_genes),
    length(overlap_seed),
    length(overlap_leaf),
    length(dr$gene_id)
  )
)



write.table(
  summary_table,
  file.path(
    data_dir,
    "rhythmicity_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)