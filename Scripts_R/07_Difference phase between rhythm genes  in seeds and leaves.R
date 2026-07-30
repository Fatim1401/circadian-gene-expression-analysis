############################################################
# Script: 07_Comparison_seeds_leaves.R
# Project: Circadian transcriptomic analysis of Pisum sativum
# Description:
# Comparison of rhythmic genes between seeds and leaves :
# - Overlap visualization (Euler diagram)
# - Phase difference analysis between organs
# - Amplitude comparison between organs
#
# Input:
# - Cycling_leaf.txt   : rhythmic genes in leaves
# - Cycling_seeds.txt  : rhythmic genes in seeds
#
# Output:
# - Euler diagram of overlap
# - Barplot of phase differences
# - Violin/density plot of amplitudes
#
# Author: Fatoumata KAMISSOKO
############################################################

## Load libraries -------------------------------------------------------------
library(tidyverse)
library(see)
library(eulerr)
library(ggpointdensity)
library(viridis)

## Load data ------------------------------------------------------------------
rhythmic_leaves <- read.table("outputs/data/Cycling_leaf.txt",  header = TRUE)
rhythmic_seeds  <- read.table("outputs/data/Cycling_seeds.txt", header = TRUE)

############################################################
# I. OVERLAP BETWEEN RHYTHMIC GENES IN SEEDS AND LEAVES
############################################################

rhythmic <- list(
  Leaves = rhythmic_leaves$CycID,
  Seeds  = rhythmic_seeds$CycID
)

plot(euler(rhythmic, shape = "ellipse"), quantities = TRUE)

############################################################
# II. PHASE DIFFERENCE BETWEEN SEEDS AND LEAVES
############################################################

# Define overlapping genes
overlapping_genes <- intersect(rhythmic_leaves$CycID, rhythmic_seeds$CycID)

# Filter phases for overlapping genes
phases_seeds_overlap  <- filter(rhythmic_seeds,  CycID %in% overlapping_genes)
phases_leaves_overlap <- filter(rhythmic_leaves, CycID %in% overlapping_genes)

phases_seeds_overlap  <- phases_seeds_overlap[,  c("CycID", "JTK_adjphase")]
phases_leaves_overlap <- phases_leaves_overlap[, c("CycID", "JTK_adjphase")]

# Rename columns
names(phases_seeds_overlap)[2]  <- "Phases_seed"
names(phases_leaves_overlap)[2] <- "Phases_leaf"

# Replace phases 24 and 26 by 0 and 2
phases_seeds_overlap$Phases_seed   <- str_replace_all(phases_seeds_overlap$Phases_seed,  "24", "0")
phases_seeds_overlap$Phases_seed   <- str_replace_all(phases_seeds_overlap$Phases_seed,  "26", "2")
phases_leaves_overlap$Phases_leaf  <- str_replace_all(phases_leaves_overlap$Phases_leaf, "24", "0")
phases_leaves_overlap$Phases_leaf  <- str_replace_all(phases_leaves_overlap$Phases_leaf, "26", "2")

# Merge phase info
phases_overlap <- merge.data.frame(phases_seeds_overlap, phases_leaves_overlap, by = "CycID")

# Calculate phase difference
phases_overlap$Phases_seed <- as.integer(phases_overlap$Phases_seed)
phases_overlap$Phases_leaf <- as.integer(phases_overlap$Phases_leaf)
phases_overlap$Phase_diff  <- phases_overlap$Phases_seed - phases_overlap$Phases_leaf

# Correct values > 12 or < -12
phases_overlap$Phase_diff[phases_overlap$Phase_diff == -22] <-  2
phases_overlap$Phase_diff[phases_overlap$Phase_diff == -20] <-  4
phases_overlap$Phase_diff[phases_overlap$Phase_diff == -18] <-  6
phases_overlap$Phase_diff[phases_overlap$Phase_diff == -16] <-  8
phases_overlap$Phase_diff[phases_overlap$Phase_diff == -14] <- 10
phases_overlap$Phase_diff[phases_overlap$Phase_diff ==  22] <- -2
phases_overlap$Phase_diff[phases_overlap$Phase_diff ==  20] <- -4
phases_overlap$Phase_diff[phases_overlap$Phase_diff ==  18] <- -6
phases_overlap$Phase_diff[phases_overlap$Phase_diff ==  16] <- -8
phases_overlap$Phase_diff[phases_overlap$Phase_diff ==  14] <- -10

# Summary statistics
cat("Median phase difference :", median(phases_overlap$Phase_diff), "\n")
cat("Mean phase difference   :", mean(phases_overlap$Phase_diff),   "\n")

count_phase_diff <- data.frame(table(phases_overlap$Phase_diff))

# Plot phase difference barplot
cairo_pdf("outputs/figures/barplot_phase_difference.pdf", width = 8, height = 6)
ggplot(count_phase_diff) +
  geom_bar(aes(x = Var1, y = Freq), stat = "identity") +
  theme_bw() +
  theme(
    axis.text         = element_text(color = "black", size = 15),
    panel.grid        = element_blank(),
    axis.ticks.length = unit(-0.15, "cm"),
    axis.text.x       = element_text(margin = margin(t = 3, unit = "mm")),
    axis.text.y       = element_text(margin = margin(r = 3, unit = "mm")),
    axis.title        = element_text(size = 15)
  ) +
  ylab("Number of rhythmic genes") +
  xlab("Phase difference (h) between seeds and leaves")
dev.off()

############################################################
# III. AMPLITUDE COMPARISON BETWEEN SEEDS AND LEAVES
############################################################

# Build amplitude table
amplitudes <- rbind(
  data.frame(
    Organ     = rep("Seeds",  length(rhythmic_seeds$CycID)),
    Amplitude = rhythmic_seeds$JTK_amplitude,
    BH.Q      = rhythmic_seeds$meta2d_BH.Q
  ),
  data.frame(
    Organ     = rep("Leaves", length(rhythmic_leaves$CycID)),
    Amplitude = rhythmic_leaves$JTK_amplitude,
    BH.Q      = rhythmic_leaves$meta2d_BH.Q
  )
)

# Calculate -log10(BH.Q)
amplitudes$log10_BH.Q <- -(log10(amplitudes$BH.Q))

# Plot amplitude by organ (violin)
cairo_pdf("outputs/figures/violin_amplitude.pdf", width = 6, height = 6)
ggplot(amplitudes, aes(x = Organ, y = Amplitude, color = Organ)) +
  geom_violinhalf(position = position_nudge(x = 0.2, y = 0)) +
  geom_jitter(alpha = 0.01, width = 0.15) +
  theme_bw() +
  theme(
    axis.text         = element_text(color = "black", size = 12),
    panel.grid        = element_blank(),
    axis.ticks.length = unit(-0.15, "cm"),
    axis.title        = element_text(size = 12)
  )
dev.off()

# Plot amplitude vs BH.Q by organ (density)
cairo_pdf("outputs/figures/density_amplitude_BHQ.pdf", width = 8, height = 5)
ggplot(amplitudes, aes(x = log10_BH.Q, y = Amplitude)) +
  geom_pointdensity() +
  scale_color_viridis() +
  geom_smooth(alpha = 0.25, color = "black", fill = "black") +
  facet_wrap(~ Organ) +
  theme_bw() +
  theme(
    axis.text         = element_text(color = "black", size = 12),
    panel.grid        = element_blank(),
    axis.ticks.length = unit(-0.15, "cm"),
    axis.text.x       = element_text(margin = margin(t = 1, unit = "mm")),
    axis.text.y       = element_text(margin = margin(r = 1, unit = "mm")),
    axis.title        = element_text(size = 12),
    strip.text        = element_text(size = 12)
  ) +
  labs(
    x = "-log10(meta2d_BH.Q)",
    y = "Amplitude of oscillation (JTK)"
  )
dev.off()



