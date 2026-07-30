############################################################
# Script: 11_Gene_Regulatory_Network_Seeds.R
# Project: Circadian transcriptomic analysis of Pisum sativum
# Description:
# Gene regulatory network inference from rhythmic transcriptome
# using dynGENIE3 algorithm — Seeds.
# - Network inference from rhythmic genes in seeds
# - Identification of rhythmic transcription factors (TFs)
# - Extraction of regulatory connections (weight threshold = 0.01)
#
# Input:
# - overlapseedsDF.txt              : rhythmic genes in seeds
# - expression.seeds.txt            : VST normalized expression (seeds)
# - TF_ITAK_v2.txt                  : transcription factor list
#
# Output:
# - link_list_dynGENIE3_seeds.txt   : list of regulatory connections
#
# Author: Fatoumata KAMISSOKO
############################################################

## Load libraries -------------------------------------------------------------
library(doRNG)
library(doParallel)
library(stringr)
library(dplyr)
library(tidyr)

## Configure dynGENIE3 --------------------------------------------------------
# Configure Rtools PATH (required for compilation on Windows)
# write('PATH="${RTOOLS45_HOME}\\usr\\bin;${PATH}"',
#       file = "~/.Renviron", append = TRUE)

# Verify Rtools installation
Sys.which("make")

# Set working directory to dynGENIE3 source folder
setwd("E:/Stage_Fatoumata/Projet-de-Stage-/data/raw/dynGENIE3-master/dynGENIE3_R_C_wrapper")

# Compile dynGENIE3 C source code (run once)
# system('R CMD SHLIB dynGENIE3.c')
# On Windows, rename .dll to .so
# file.rename("dynGENIE3.dll", "dynGENIE3.so")

# Load dynGENIE3
source("dynGENIE3.R")
dyn.load("dynGENIE3.so")

## Load data ------------------------------------------------------------------
cycling_TOT_seed <- read.table("outputs/data/overlapseedsDF.txt",    header = TRUE)
vst100           <- read.table("outputs/data/expression.seeds.txt",  header = TRUE)
TF               <- read.table("E:/Stage_Fatoumata/TF_ITAK_v2.txt", header = TRUE)

############################################################
# I. PREPARE EXPRESSION DATA
############################################################

# Filter expression data for rhythmic genes
data_rhythmic <- filter(vst100, ID %in% cycling_TOT_seed$gene_id)

# Define row names and add time points row
rownames(data_rhythmic)     <- data_rhythmic[, 1]
data_rhythmic               <- rbind(names(data_rhythmic), data_rhythmic)
data_rhythmic               <- data_rhythmic[, -1]
row.names(data_rhythmic)[1] <- "time.points"

# Assign numeric time points (0 replaced by 0.001 — dynGENIE3 requirement)
data_rhythmic[1, c("S1",  "S2",  "S3")]  <- 0.001  # T0
data_rhythmic[1, c("S5",  "S6",  "S7")]  <- 4      # T4
data_rhythmic[1, c("S9",  "S10", "S11")] <- 8      # T8
data_rhythmic[1, c("S13", "S14", "S15")] <- 12     # T12
data_rhythmic[1, c("S17", "S18", "S19")] <- 16     # T16
data_rhythmic[1, c("S21", "S22", "S23")] <- 20     # T20
data_rhythmic[1, c("S25", "S26", "S27")] <- 24     # T24
data_rhythmic[1, c("S29", "S30", "S31")] <- 28     # T28
data_rhythmic[1, c("S33", "S34", "S35")] <- 32     # T32
data_rhythmic[1, c("S37", "S38", "S39")] <- 36     # T36
data_rhythmic[1, c("S41", "S42", "S43")] <- 40     # T40
data_rhythmic[1, c("S45", "S46", "S47")] <- 44     # T44
data_rhythmic[1, c("S49", "S50", "S51")] <- 48     # T48

data_rhythmic <- data_rhythmic %>% mutate_if(is.character, as.numeric)

############################################################
# II. SEPARATE DATA BY REPLICATE
############################################################

rep1 <- data_rhythmic[, c("S1","S5","S9","S13","S17","S21","S25","S29","S33","S37","S41","S45","S49")]
rep2 <- data_rhythmic[, c("S2","S6","S10","S14","S18","S22","S26","S30","S34","S38","S42","S46","S50")]
rep3 <- data_rhythmic[, c("S3","S7","S11","S15","S19","S23","S27","S31","S35","S39","S43","S47","S51")]

cat("Dimensions rep1 :", dim(rep1), "\n")

# Check for duplicated values (dynGENIE3 requirement)
rep1_dup <- apply(rep1, 1, FUN = function(x) ifelse(max(table(x)) > (ncol(rep1) - 1), TRUE, FALSE))
rep2_dup <- apply(rep2, 1, FUN = function(x) ifelse(max(table(x)) > (ncol(rep2) - 1), TRUE, FALSE))
rep3_dup <- apply(rep3, 1, FUN = function(x) ifelse(max(table(x)) > (ncol(rep3) - 1), TRUE, FALSE))

cat("Genes to remove — rep1 :", sum(rep1_dup), "\n")
cat("Genes to remove — rep2 :", sum(rep2_dup), "\n")
cat("Genes to remove — rep3 :", sum(rep3_dup), "\n")

############################################################
# III. DEFINE REGULATORS (RHYTHMIC TFs)
############################################################

rhythmic_TF <- filter(TF, ID %in% cycling_TOT_seed$gene_id)
input.genes <- unique(rhythmic_TF$ID)

cat("Number of rhythmic TFs :", length(input.genes), "\n")

############################################################
# IV. RUN dynGENIE3
############################################################

# Convert to matrix
rep1 <- as.matrix(rep1)
rep2 <- as.matrix(rep2)
rep3 <- as.matrix(rep3)

# Reformat for dynGENIE3
time.points <- list(rep1[1, ], rep2[1, ], rep3[1, ])
TS.data     <- list(
  rep1[2:nrow(rep1), ],
  rep2[2:nrow(rep2), ],
  rep3[2:nrow(rep3), ]
)

# Run network inference
res        <- dynGENIE3(TS.data, time.points, regulators = input.genes)
res.weight <- as.data.frame(res$weight.matrix)

# Extract regulatory connections (weight threshold = 0.01)
link.list <- get.link.list(res$weight.matrix, threshold = 0.01)

cat("Number of regulatory connections :", nrow(link.list), "\n")

############################################################
# V. EXPORT RESULTS
############################################################

write.table(link.list,
            "outputs/data/link_list_dynGENIE3_seeds.txt",
            row.names = FALSE,
            quote     = FALSE,
            sep       = "\t")