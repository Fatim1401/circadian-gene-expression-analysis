############################################################
# Script: 12_Network_GO_Filtering.R
# Project: Circadian transcriptomic analysis of Pisum sativum
# Description:
# Filtering of gene regulatory networks based on GO annotation.
# - Filter regulatory connections by weight threshold
# - Extract sub-networks for abiotic stress response (GO:0009628)
#   in leaves and seeds
# - Count regulators and target genes per organ
#
# Input:
# - link_list_dynGENIE3.txt        : regulatory network (leaves)
# - link_list_dynGENIE3_seed.txt   : regulatory network (seeds)
# - GO_description.csv             : GO term descriptions
# - Gene_annotation_v2.csv         : GO annotations per gene
#
# Output:
# - GO_network_abiotic_Target_L.txt : abiotic stress network (leaves)
# - GO_network_abiotic_Target_S.txt : abiotic stress network (seeds)
#
# Author: Fatoumata KAMISSOKO
############################################################

## Load libraries -------------------------------------------------------------
library(dplyr)

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
link.list_leaf30 <- link.list_leaf %>%
  filter(weight >= 0.03)

# Extract GO annotations for abiotic stress genes
Go_annot_abiotic_L <- GO_annotation_v2 %>%
  filter(GO %in% Go_abiotic_L)

# Filter network for abiotic stress target genes
GO_network_abiotic_Target_L <- link.list_leaf30 %>%
  filter(target.gene %in% Go_annot_abiotic_L$GID)

# Summary
cat("Leaves — Unique regulators :", length(unique(GO_network_abiotic_Target_L$regulatory.gene)), "\n")
cat("Leaves — Unique targets    :", length(unique(GO_network_abiotic_Target_L$target.gene)), "\n")

# Count connections per regulator
count.regulator_leaf <- data.frame(table(GO_network_abiotic_Target_L$regulatory.gene))
colnames(count.regulator_leaf) <- c("regulatory.gene", "n_targets")
count.regulator_leaf <- count.regulator_leaf[order(count.regulator_leaf$n_targets,
                                                   decreasing = TRUE), ]

############################################################
# II. NETWORK FILTERING — SEEDS
############################################################

# Define GO term of interest
Go_abiotic_S <- "GO:0009628"

# Extract GO annotations for abiotic stress genes
Go_annot_abiotic_S <- GO_annotation_v2 %>%
  filter(GO %in% Go_abiotic_S)

# Filter network for abiotic stress target genes
GO_network_abiotic_Target_S <- link.list_seed %>%
  filter(target.gene %in% Go_annot_abiotic_S$GID)

# Summary
cat("Seeds — Unique regulators :", length(unique(GO_network_abiotic_Target_S$regulatory.gene)), "\n")
cat("Seeds — Unique targets    :", length(unique(GO_network_abiotic_Target_S$target.gene)), "\n")

# Count connections per regulator
count.regulator_seed <- data.frame(table(GO_network_abiotic_Target_S$regulatory.gene))
colnames(count.regulator_seed) <- c("regulatory.gene", "n_targets")
count.regulator_seed <- count.regulator_seed[order(count.regulator_seed$n_targets,
                                                   decreasing = TRUE), ]

############################################################
# III. EXPORT RESULTS
############################################################

write.table(GO_network_abiotic_Target_L,
            "outputs/data/GO_network_abiotic_Target_L.txt",
            row.names = FALSE,
            quote     = FALSE,
            sep       = "\t")

write.table(GO_network_abiotic_Target_S,
            "outputs/data/GO_network_abiotic_Target_S.txt",
            row.names = FALSE,
            quote     = FALSE,
            sep       = "\t")