############################################################

# Script: 06A_PCA_analysis.R
# Project: Circadian transcriptomic analysis under abiotic stress

# Description:
# Principal Component Analysis (PCA) of normalized expression
# data for seeds and leaves. Visualization includes scree plots
# and PCA projections (PC1–PC2 and PC3–PC4) colored by time.
#

# Input:
# - Expression matrices (seeds and leaves)
# - Metadata file
#

# Output:
# - PCA plots (PDF)
#

# Author: Fatoumata KAMISSOKO 

############################################################

## Load libraries -------------------------------------------------------------

library(ade4)
library(factoextra)
library(tidyverse)
library(viridis)
library(ggalt)

## Define paths ---------------------------------------------------------------

metadata_file        <- "data/raw/Metadata.txt"
expression_seeds_file  <- "outputs/data/expression.seeds.txt"
expression_leaves_file <- "outputs/data/expression.leaves.txt"
output_dir          <- "outputs/figures"

## Load data ------------------------------------------------------------------

metadata <- read.table(metadata_file, header = TRUE, sep = "\t")
expression.seeds  <- read.table(expression_seeds_file, header = TRUE)
expression.leaves <- read.table(expression_leaves_file, header = TRUE)

############################################################

# I. PCA — SEEDS

############################################################

row.names(expression.seeds) <- expression.seeds$ID
expr_seeds <- as.data.frame(t(expression.seeds[,-1]))

expr_seeds$SampleName <- row.names(expr_seeds)
expr_seeds <- merge(metadata, expr_seeds, by = "SampleName")
row.names(expr_seeds) <- expr_seeds$SampleName

res.pca <- dudi.pca(expr_seeds[,-c(1:4)], nf = 5, scannf = FALSE)
eig.val <- get_eigenvalue(res.pca)

# Scree plot

p_scree <- fviz_screeplot(res.pca, addlabels = TRUE, ncp = 10) +
  theme_bw(base_size = 18)

ggsave(file.path(output_dir, "Seeds_Screeplot.pdf"), p_scree, width = 7, height = 6)

# Axis labels

PC1_title <- paste0("PC1 (", round(eig.val[1,2], 2), "%)")
PC2_title <- paste0("PC2 (", round(eig.val[2,2], 2), "%)")
PC3_title <- paste0("PC3 (", round(eig.val[3,2], 2), "%)")
PC4_title <- paste0("PC4 (", round(eig.val[4,2], 2), "%)")

expr_seeds$Time <- factor(expr_seeds$Time,
                          levels = paste0("ZT", seq(0,48,4)))

df_pca <- cbind(res.pca$li, Time = expr_seeds$Time)
time_n <- length(unique(expr_seeds$Time))

p_base <- ggplot() + coord_fixed() +
  scale_color_manual(values = viridis(time_n)) +
  theme_bw()

# PC1-PC2

p12 <- p_base +
  geom_point(data = df_pca, aes(Axis1, Axis2, color = Time), size = 3) +
  geom_encircle(data = df_pca, aes(Axis1, Axis2, group = Time, color = Time), alpha = 0.25) +
  labs(x = PC1_title, y = PC2_title)

ggsave(file.path(output_dir, "Seeds_PCA_PC1_PC2.pdf"), p12, width = 7, height = 6)

# PC3-PC4

p34 <- p_base +
  geom_point(data = df_pca, aes(Axis3, Axis4, color = Time), size = 3) +
  geom_encircle(data = df_pca, aes(Axis3, Axis4, group = Time, color = Time), alpha = 0.25) +
  labs(x = PC3_title, y = PC4_title)

ggsave(file.path(output_dir, "Seeds_PCA_PC3_PC4.pdf"), p34, width = 7, height = 6)

############################################################

# II. PCA — LEAVES

############################################################

row.names(expression.leaves) <- expression.leaves$ID
expr_leaves <- as.data.frame(t(expression.leaves[,-1]))

expr_leaves$SampleName <- row.names(expr_leaves)
expr_leaves <- merge(metadata, expr_leaves, by = "SampleName")
row.names(expr_leaves) <- expr_leaves$SampleName

res.pca <- dudi.pca(expr_leaves[,-c(1:4)], nf = 5, scannf = FALSE)
eig.val <- get_eigenvalue(res.pca)

# Scree plot

p_scree <- fviz_screeplot(res.pca, addlabels = TRUE, ncp = 10) +
  theme_bw(base_size = 18)

ggsave(file.path(output_dir, "Leaves_Screeplot.pdf"), p_scree, width = 7, height = 6)

# Axis labels

PC1_title <- paste0("PC1 (", round(eig.val[1,2], 2), "%)")
PC2_title <- paste0("PC2 (", round(eig.val[2,2], 2), "%)")
PC3_title <- paste0("PC3 (", round(eig.val[3,2], 2), "%)")
PC4_title <- paste0("PC4 (", round(eig.val[4,2], 2), "%)")

expr_leaves$Time <- factor(expr_leaves$Time,
                           levels = paste0("ZT", seq(0,48,4)))

df_pcaL <- cbind(res.pca$li, Time = expr_leaves$Time)
time_n <- length(unique(expr_leaves$Time))

# PC1-PC2

p12 <- p_base +
  geom_point(data = df_pcaL, aes(Axis1, Axis2, color = Time), size = 3) +
  geom_encircle(data = df_pcaL, aes(Axis1, Axis2, group = Time, color = Time), alpha = 0.25) +
  labs(x = PC1_title, y = PC2_title)

ggsave(file.path(output_dir, "Leaves_PCA_PC1_PC2.pdf"), p12, width = 7, height = 6)

# PC3-PC4

p34 <- p_base +
  geom_point(data = df_pcaL, aes(Axis3, Axis4, color = Time), size = 3) +
  geom_encircle(data = df_pcaL, aes(Axis3, Axis4, group = Time, color = Time), alpha = 0.25) +
  labs(x = PC3_title, y = PC4_title)

ggsave(file.path(output_dir, "Leaves_PCA_PC3_PC4.pdf"), p34, width = 7, height = 6)
