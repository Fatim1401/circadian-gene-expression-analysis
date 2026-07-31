############################################################
# Script: 11_Network_Phase_Analysis.R
# Project: Circadian transcriptomic analysis of Pisum sativum
# Description:
# Phase analysis of regulators and target genes in abiotic
# stress gene regulatory networks for leaves and seeds.
# - Filter networks by weight threshold (>= 0.025)
# - Extract sub-networks where both regulators and targets
#   are annotated with abiotic stress GO term (GO:0009628)
# - Phase distribution of TFs and target genes
# - Top 15 TF barplots and overlap between organs
#
# Input:
# - link_list_dynGENIE3.txt         : regulatory network (leaves)
# - link_list_dynGENIE3_seed.txt    : regulatory network (seeds)
# - GO_description.csv              : GO term descriptions
# - Gene_annotation_v2.csv          : GO annotations per gene
# - Leaf_Metacycle_results.txt      : phase info (leaves)
# - Seed_Metacycle_results.txt      : phase info (seeds)
# - TOP_15_TF_S1.xlsx               : Top 15 TFs seeds
# - Top_15_TF_1.xlsx                : Top 15 TFs leaves
#
# Output:
# - filter_reg_tar.txt              : filtered network (leaves)
# - filter_reg_tar_s.txt            : filtered network (seeds)
# - count.regulator1.txt/.xlsx      : regulator counts (leaves)
# - count.regulator_s.txt/.xlsx     : regulator counts (seeds)
# - Phase_TF_leaf.pdf               : phase barplot TFs leaves
# - Phase_Target_leaf.pdf           : phase barplot targets leaves
# - Phase_TF_seed.pdf               : phase barplot TFs seeds
# - Phase_Target_seed.pdf           : phase barplot targets seeds
# - top15_graine.pdf                : Top 15 TFs seeds barplot
# - top15_feuilles.pdf              : Top 15 TFs leaves barplot
# - Venn_Graine_feuilles_top15.pdf  : Venn diagram Top 15 TFs
#
# Author: Fatoumata KAMISSOKO
############################################################

## Load libraries -------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(openxlsx)
library(readxl)
library(VennDiagram)
library(grid)

## Load data ------------------------------------------------------------------
link.list_leaf   <- read.table("outputs/data/link_list_dynGENIE3.txt",
                               header = TRUE)
link.list_seed   <- read.table("outputs/data/link_list_dynGENIE3_seed.txt",
                               header = TRUE)
GO_description   <- read.csv("outputs/data/GO_description.csv",
                             header = TRUE)
GO_annotation_v2 <- read.csv2("outputs/data/Gene_annotation_v2.csv",
                              header = TRUE)

############################################################
# I. NETWORK FILTERING — LEAVES
############################################################

# Define GO term of interest
Go_abiotic_L <- "GO:0009628"

# Filter network by weight threshold
link.list_leaf_0.025 <- link.list_leaf %>%
  filter(weight >= 0.025)

# Extract GO annotations for abiotic stress
Go_annot_abiotic_L <- GO_annotation_v2 %>%
  filter(GO %in% Go_abiotic_L)

# Filter target genes annotated with abiotic stress GO
GO_network_abiotic_Target_L <- link.list_leaf_0.025 %>%
  filter(target.gene %in% Go_annot_abiotic_L$GID)

# Filter regulatory genes also annotated with abiotic stress GO
filter_reg_tar <- GO_network_abiotic_Target_L %>%
  filter(regulatory.gene %in% Go_annot_abiotic_L$GID)

# Summary
cat("Leaves — Unique regulators :", length(unique(filter_reg_tar$regulatory.gene)), "\n")
cat("Leaves — Unique targets    :", length(unique(filter_reg_tar$target.gene)), "\n")

# Count connections per regulator
count.regulator <- data.frame(table(filter_reg_tar$regulatory.gene))
colnames(count.regulator) <- c("regulatory.gene", "n_targets")
count.regulator <- count.regulator[order(count.regulator$n_targets, decreasing = TRUE), ]

## Export data ----------------------------------------------------------------
write.table(filter_reg_tar,
            "outputs/data/filter_reg_tar.txt",
            sep       = "\t",
            row.names = FALSE,
            quote     = FALSE)

write.table(count.regulator,
            "outputs/data/count.regulator1.txt",
            sep       = "\t",
            row.names = FALSE,
            quote     = FALSE)

write.xlsx(count.regulator,
           "outputs/data/count.regulator1.xlsx",
           rowNames = FALSE)

############################################################
# II. PHASE ANALYSIS — LEAVES
############################################################

phase_table_leaf <- read.table("outputs/data/Leaf_Metacycle_results.txt",
                               header = TRUE)

# Phase correction function
correct_phase <- function(df, phase_col = "JTK_adjphase") {
  df %>% mutate(!!phase_col := case_when(
    !!sym(phase_col) == 24 ~ 0,
    !!sym(phase_col) == 26 ~ 2,
    TRUE ~ !!sym(phase_col)
  ))
}

# TF phase distribution
merged_data_TF <- filter_reg_tar %>%
  left_join(phase_table_leaf, by = c("regulatory.gene" = "CycID")) %>%
  correct_phase()

TF_unique <- merged_data_TF %>%
  distinct(regulatory.gene, JTK_adjphase)

count.regulator_phase_TF <- data.frame(table(TF_unique$JTK_adjphase))
colnames(count.regulator_phase_TF) <- c("phase", "n")

cairo_pdf("outputs/figures/Phase_TF_leaf.pdf", width = 9, height = 9)
ggplot(count.regulator_phase_TF, aes(x = phase, y = n)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.3, size = 5) +
  theme_bw() +
  theme(
    axis.text         = element_text(color = "black", size = 15),
    panel.grid        = element_blank(),
    axis.ticks.length = unit(-0.15, "cm"),
    axis.text.x       = element_text(margin = margin(3, 0, 0, 0, unit = "mm")),
    axis.text.y       = element_text(margin = margin(0, 3, 0, 0, unit = "mm")),
    axis.title        = element_text(size = 15)
  ) +
  labs(
    subtitle = paste0("Total: ", sum(count.regulator_phase_TF$n, na.rm = TRUE),
                      " Transcription factors"),
    x = "Phase (ZT)",
    y = "Number of regulatory genes"
  )
dev.off()

# Target phase distribution
merged_data_target <- filter_reg_tar %>%
  left_join(phase_table_leaf, by = c("target.gene" = "CycID")) %>%
  correct_phase()

Target_unique <- merged_data_target %>%
  distinct(target.gene, JTK_adjphase)

count.regulator_phase_target <- data.frame(table(Target_unique$JTK_adjphase))
colnames(count.regulator_phase_target) <- c("phase", "n")

cairo_pdf("outputs/figures/Phase_Target_leaf.pdf", width = 9, height = 9)
ggplot(count.regulator_phase_target, aes(x = phase, y = n)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.3, size = 5) +
  theme_bw() +
  theme(
    axis.text         = element_text(color = "black", size = 15),
    panel.grid        = element_blank(),
    axis.ticks.length = unit(-0.15, "cm"),
    axis.text.x       = element_text(margin = margin(3, 0, 0, 0, unit = "mm")),
    axis.text.y       = element_text(margin = margin(0, 3, 0, 0, unit = "mm")),
    axis.title        = element_text(size = 15)
  ) +
  labs(
    subtitle = paste0("Total: ", sum(count.regulator_phase_target$n, na.rm = TRUE),
                      " Target genes"),
    x = "Phase (ZT)",
    y = "Number of target genes"
  )
dev.off()

############################################################
# III. NETWORK FILTERING — SEEDS
############################################################

# Define GO term of interest
Go_abiotic_S <- "GO:0009628"

# Filter network by weight threshold
link.list_seed_0.025 <- link.list_seed %>%
  filter(weight >= 0.025)

# Extract GO annotations for abiotic stress
Go_annot_abiotic_S <- GO_annotation_v2 %>%
  filter(GO %in% Go_abiotic_S)

# Filter target genes annotated with abiotic stress GO
GO_network_abiotic_Target_S <- link.list_seed_0.025 %>%
  filter(target.gene %in% Go_annot_abiotic_S$GID)

# Filter regulatory genes also annotated with abiotic stress GO
filter_reg_tar_s <- GO_network_abiotic_Target_S %>%
  filter(regulatory.gene %in% Go_annot_abiotic_S$GID)

# Summary
cat("Seeds — Unique regulators :", length(unique(filter_reg_tar_s$regulatory.gene)), "\n")
cat("Seeds — Unique targets    :", length(unique(filter_reg_tar_s$target.gene)), "\n")

# Count connections per regulator
count.regulator_s <- data.frame(table(filter_reg_tar_s$regulatory.gene))
colnames(count.regulator_s) <- c("regulatory.gene", "n_targets")
count.regulator_s <- count.regulator_s[order(count.regulator_s$n_targets, decreasing = TRUE), ]

## Export data ----------------------------------------------------------------
write.table(filter_reg_tar_s,
            "outputs/data/filter_reg_tar_s.txt",
            sep       = "\t",
            row.names = FALSE,
            quote     = FALSE)

write.table(count.regulator_s,
            "outputs/data/count.regulator_s.txt",
            sep       = "\t",
            row.names = FALSE,
            quote     = FALSE)

write.xlsx(count.regulator_s,
           "outputs/data/count.regulator_s.xlsx",
           rowNames = FALSE)

############################################################
# IV. PHASE ANALYSIS — SEEDS
############################################################

phase_table_seed <- read.table("outputs/data/Seed_Metacycle_results.txt",
                               header = TRUE)

# TF phase distribution
merged_data_TF_s <- filter_reg_tar_s %>%
  left_join(phase_table_seed, by = c("regulatory.gene" = "CycID")) %>%
  correct_phase()

TF_unique_s <- merged_data_TF_s %>%
  distinct(regulatory.gene, JTK_adjphase)

count.regulator_phase_TF_s <- data.frame(table(TF_unique_s$JTK_adjphase))
colnames(count.regulator_phase_TF_s) <- c("phase", "n")

cairo_pdf("outputs/figures/Phase_TF_seed.pdf", width = 9, height = 9)
ggplot(count.regulator_phase_TF_s, aes(x = phase, y = n)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.3, size = 5) +
  theme_bw() +
  theme(
    axis.text         = element_text(color = "black", size = 15),
    panel.grid        = element_blank(),
    axis.ticks.length = unit(-0.15, "cm"),
    axis.text.x       = element_text(margin = margin(3, 0, 0, 0, unit = "mm")),
    axis.text.y       = element_text(margin = margin(0, 3, 0, 0, unit = "mm")),
    axis.title        = element_text(size = 15)
  ) +
  labs(
    subtitle = paste0("Total: ", sum(count.regulator_phase_TF_s$n, na.rm = TRUE),
                      " Transcription factors"),
    x = "Phase (ZT)",
    y = "Number of regulatory genes"
  )
dev.off()

# Target phase distribution
merged_data_target_s <- filter_reg_tar_s %>%
  left_join(phase_table_seed, by = c("target.gene" = "CycID")) %>%
  correct_phase()

Target_unique_s <- merged_data_target_s %>%
  distinct(target.gene, JTK_adjphase)

count.regulator_phase_target_s <- data.frame(table(Target_unique_s$JTK_adjphase))
colnames(count.regulator_phase_target_s) <- c("phase", "n")

cairo_pdf("outputs/figures/Phase_Target_seed.pdf", width = 9, height = 9)
ggplot(count.regulator_phase_target_s, aes(x = phase, y = n)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.3, size = 5) +
  theme_bw() +
  theme(
    axis.text         = element_text(color = "black", size = 15),
    panel.grid        = element_blank(),
    axis.ticks.length = unit(-0.15, "cm"),
    axis.text.x       = element_text(margin = margin(3, 0, 0, 0, unit = "mm")),
    axis.text.y       = element_text(margin = margin(0, 3, 0, 0, unit = "mm")),
    axis.title        = element_text(size = 15)
  ) +
  labs(
    subtitle = paste0("Total: ", sum(count.regulator_phase_target_s$n, na.rm = TRUE),
                      " Target genes"),
    x = "Phase (ZT)",
    y = "Number of target genes"
  )
dev.off()

############################################################
# V. TOP 15 TFs BARPLOTS AND OVERLAP
############################################################

# Load Top 15 TF tables
Top15_seed <- read_excel("outputs/data/TOP_15_TF_S1.xlsx")
Top15_leaf <- read_excel("outputs/data/Top_15_TF_1.xlsx")

# Add organ label and select columns
Top15_seed <- Top15_seed %>%
  mutate(condition = "Seeds") %>%
  select(IDPS_V2, IDPS_V1, DESCRIPTION_PS, Freq, condition)

Top15_leaf <- Top15_leaf %>%
  mutate(condition = "Leaves") %>%
  select(IDPS_V2, IDPS_V1, DESCRIPTION_PS, Freq, condition)

# Seeds barplot
cairo_pdf("outputs/figures/top15_seeds.pdf", width = 6, height = 7)
ggplot(Top15_seed, aes(x = reorder(DESCRIPTION_PS, Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "#e31a1c", width = 0.7) +
  coord_flip() +
  labs(
    title = "Top 15 Transcription factors — Seeds",
    x     = "Transcription factors",
    y     = "Number of regulated genes"
  ) +
  theme_classic() +
  theme(
    plot.title  = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10)
  )
dev.off()

# Leaves barplot
cairo_pdf("outputs/figures/top15_leaves.pdf", width = 6, height = 7)
ggplot(Top15_leaf, aes(x = reorder(DESCRIPTION_PS, Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "#1f78b4", width = 0.7) +
  coord_flip() +
  labs(
    title = "Top 15 Transcription factors — Leaves",
    x     = "Transcription factors",
    y     = "Number of regulated genes"
  ) +
  theme_classic() +
  theme(
    plot.title  = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10)
  )
dev.off()

# Venn diagram — overlap between Top 15 TFs
overlap_leaf_seed <- intersect(Top15_seed$DESCRIPTION_PS, Top15_leaf$DESCRIPTION_PS)

cat("Common TFs between seeds and leaves :", length(overlap_leaf_seed), "\n")
print(overlap_leaf_seed)

cairo_pdf("outputs/figures/Venn_seeds_leaves_top15.pdf", width = 6, height = 6)
draw.pairwise.venn(
  area1      = length(Top15_leaf$DESCRIPTION_PS),
  area2      = length(Top15_seed$DESCRIPTION_PS),
  cross.area = length(overlap_leaf_seed),
  category   = c("Leaves", "Seeds"),
  fill       = c("#1f78b4", "#e31a1c"),
  alpha      = 0.7,
  cex        = 1.3,
  cat.cex    = 1.2,
  lty        = "blank"
)
grid.text(
  paste(overlap_leaf_seed, collapse = "\n"),
  x  = 0.5,
  y  = 0.5,
  gp = gpar(fontsize = 8)
)
dev.off()