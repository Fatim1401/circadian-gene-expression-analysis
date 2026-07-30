############################################################

# Script: 09_Individual_gene_expression_plot.R
# Project: Circadian transcriptomic GENE expression analysis

# Description:

# Visualization of individual gene expression profiles across
# circadian time in seeds and leaves, including:
# - Mean ± SD expression
# - Cosine fitting
# - Phase comparison between organs
#

# Input:
# - Expression matrices (seeds & leaves)
# - Metadata
# - MetaCycle results (phase information)
#

# Output:
# - Individual gene expression plots (PDF)

#

# Author: Fatoumata KAMISSOKO 

############################################################

## Load libraries -------------------------------------------------------------

library(tidyverse)
library(reshape2)
library(ggh4x)

## Load data ------------------------------------------------------------------

expression.seeds <- read.table("outputs/data/expression.seeds.txt", header = TRUE)
expression.leaves <- read.table("outputs/data/expression.leaves.txt", header = TRUE)
metadata <- read.table("data/raw/Metadata.txt", header = TRUE)

## Select gene ---------------------------------------------------------------

Gene.ID <- "Psat0s580g0080"

############################################################

# I. Prepare expression data

############################################################

expression.seeds <- filter(expression.seeds, ID %in% Gene.ID)
expression.leaves <- filter(expression.leaves, ID %in% Gene.ID)

expression.seeds <- melt(expression.seeds, id.vars = "ID")
expression.leaves <- melt(expression.leaves, id.vars = "ID")

expression.selected <- rbind(expression.seeds, expression.leaves)

names(metadata)[1] <- "variable"
expression.selected <- left_join(metadata, expression.selected, by = "variable")

############################################################

# II. Compute mean and SD

############################################################

data.mean <- expression.selected %>%
  group_by(Organ, Time) %>%
  summarise(Mean = mean(value), SD = sd(value), .groups = "drop")

data.mean$Time <- str_replace_all(data.mean$Time, "ZT","")
data.mean$Time <- as.numeric(data.mean$Time)
data.mean$Organ <- as.factor(data.mean$Organ)

############################################################

# III. Basic expression plot

############################################################

colors <- c("#D08D20","#424DB3")

cairo_pdf("outputs/figures/gene_expression_basic.pdf", width = 6, height = 6)

ggplot(data.mean, aes(x = Time, y = Mean, color = Organ)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_bw() +
  labs(title = Gene.ID,
       x = "Time (ZT)",
       y = "Expression (VST)")

dev.off()

############################################################

# IV. Phase comparison (MetaCycle)

############################################################

metaSeeds <- read.table("outputs/data/Seed_Metacycle_results.txt", header = TRUE)
metaLeaves <- read.table("outputs/data/Leaf_Metacycle_results.txt", header = TRUE)

phaseS <- metaSeeds[metaSeeds$CycID == Gene.ID, "JTK_adjphase"]
phaseL <- metaLeaves[metaLeaves$CycID == Gene.ID, "JTK_adjphase"]

phase.df <- data.frame(
  Organ = c("S","L"),
  Phase = as.numeric(c(phaseS, phaseL))
)

phase.diff <- phase.df$Phase[1] - phase.df$Phase[2]

############################################################

# V. Cosine fit plot

############################################################

cairo_pdf("outputs/figures/gene_expression_cosine_fit.pdf", width = 6, height = 6)

ggplot(data.mean, aes(x = Time, y = Mean, color = Organ)) +
  geom_point() +
  geom_line() +
  geom_smooth(method="lm", se=FALSE,
              formula = y ~ sin(2*pi*x/24) + cos(2*pi*x/24)) +
  geom_vline(data = phase.df, aes(xintercept = Phase, color = Organ),
             linetype = "dashed") +
  theme_bw() +
  labs(title = paste(Gene.ID, "- Phase diff:", phase.diff))

dev.off()

############################################################

# VI. Scaled comparison (Seeds vs Leaves)

############################################################

expression.seeds$value_scaled <- scale(expression.seeds$value)
expression.leaves$value_scaled <- scale(expression.leaves$value)

expr.scaled <- rbind(expression.seeds, expression.leaves)
expr.scaled <- left_join(metadata, expr.scaled, by = "variable")

data.scaled <- expr.scaled %>%
  group_by(Organ, Time) %>%
  summarise(Mean = mean(value_scaled), SD = sd(value_scaled), .groups = "drop")

data.scaled$Time <- as.numeric(str_replace_all(data.scaled$Time, "ZT",""))

cairo_pdf("outputs/figures/gene_expression_scaled.pdf", width = 6, height = 6)

ggplot(data.scaled, aes(x = Time, y = Mean, color = Organ)) +
  geom_point() +
  geom_line() +
  geom_smooth(method="lm", se=FALSE,
              formula = y ~ sin(2*pi*x/24) + cos(2*pi*x/24)) +
  theme_bw() +
  labs(title = paste(Gene.ID, "(scaled)"))

dev.off()
