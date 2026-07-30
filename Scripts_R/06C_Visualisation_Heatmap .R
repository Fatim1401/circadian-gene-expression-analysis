############################################################

# Script: 08_Heatmap_rhythmic_genes.R
# Project: Circadian transcriptomic analysis

# Description:

# Heatmap visualization of rhythmic genes identified from
# MetaCycle + LimoRhyde overlap, ordered by circadian phase.

#
# Input:
# - Expression matrices
# - Overlap gene lists with phase

# Output:
# - Heatmaps for seeds and leaves

# Author: Fatoumata KAMISSOKO 

############################################################

## Load libraries -------------------------------------------------------------

library(dplyr)
library(pheatmap)
library(ggplotify)

## Load data ------------------------------------------------------------------

expression.seeds <- read.table("outputs/data/expression.seeds.txt", header = TRUE)
expression.leaves <- read.table("outputs/data/expression.leaves.txt", header = TRUE)

overlapSeedsDF <- read.table("outputs/data/overlapSeedsDF.csv", header = TRUE)
overlapLeafDF  <- read.table("outputs/data/overlapLeafDF.csv", header = TRUE)

############################################################

# I. HEATMAP SEEDS

############################################################

Rhythmic.seeds <- overlapSeedsDF[,c("gene_id","JTK_adjphase")]
names(Rhythmic.seeds) <- c("ID","Phase")

expr <- expression.seeds %>%
  filter(ID %in% Rhythmic.seeds$ID) %>%
  left_join(Rhythmic.seeds, by="ID")

expr$Phase <- as.numeric(expr$Phase) %% 24
expr <- expr %>% arrange(Phase)

annotation_row <- expr %>% select(Phase)
rownames(annotation_row) <- expr$ID

expr_mat <- expr %>% select(-ID,-Phase) %>% as.matrix()
rownames(expr_mat) <- expr$ID

p <- pheatmap(expr_mat,
              scale="row",
              show_rownames=FALSE,
              cluster_cols=FALSE,
              cluster_rows=FALSE)

ggsave("outputs/figures/heatmap_seeds.pdf", as.ggplot(p))

############################################################

# II. HEATMAP LEAVES

############################################################

Rhythmic.leaves <- overlapLeafDF[,c("gene_id","JTK_adjphase")]
names(Rhythmic.leaves) <- c("ID","Phase")

expr <- expression.leaves %>%
  filter(ID %in% Rhythmic.leaves$ID) %>%
  left_join(Rhythmic.leaves, by="ID")

expr$Phase <- as.numeric(expr$Phase) %% 24
expr <- expr %>% arrange(Phase)

annotation_row <- expr %>% select(Phase)
rownames(annotation_row) <- expr$ID

expr_mat <- expr %>% select(-ID,-Phase) %>% as.matrix()
rownames(expr_mat) <- expr$ID

p <- pheatmap(expr_mat,
              scale="row",
              show_rownames=FALSE,
              cluster_cols=FALSE,
              cluster_rows=FALSE)

ggsave("outputs/figures/heatmap_leaves.pdf", as.ggplot(p))
