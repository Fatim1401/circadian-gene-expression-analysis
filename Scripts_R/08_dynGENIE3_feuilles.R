############################################################
# Script: 08_Gene_Regulatory_Network.R
# Project: Circadian transcriptomic analysis of Pisum sativum
# Description:
# Gene regulatory network inference from rhythmic transcriptome
# using dynGENIE3 algorithm.
# - Network inference from rhythmic genes in leaves
# - Identification of rhythmic transcription factors (TFs)
# - Extraction of regulatory connections (weight threshold = 0.01)
#
# Input:
# - overlapLeafDF.txt         : rhythmic genes in leaves
# - expression.leaves.txt     : VST normalized expression (leaves)
# - TF_ITAK_v2.txt            : transcription factor list
#
# Output:
# - link_list_dynGENIE3.txt   : list of regulatory connections
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
cycling_TOT_leaf <- read.table("outputs/data/overlapLeafDF.txt",    header = TRUE)
vst100           <- read.table("outputs/data/expression.leaves.txt", header = TRUE)
TF               <- read.table("E:/Stage_Fatoumata/TF_ITAK_v2.txt", header = TRUE)

############################################################
# I. PREPARE EXPRESSION DATA
############################################################

# Filter expression data for rhythmic genes
data_rhythmic <- filter(vst100, ID %in% cycling_TOT_leaf$AGI)

# Define row names and add time points row
rownames(data_rhythmic)      <- data_rhythmic[, 1]
data_rhythmic                <- rbind(names(data_rhythmic), data_rhythmic)
data_rhythmic                <- data_rhythmic[, -1]
row.names(data_rhythmic)[1]  <- "time.points"

# Assign numeric time points (0 replaced by 0.001 — dynGENIE3 requirement)
data_rhythmic[1, c("L1",  "L2",  "L3")]  <- 0.001  # T0
data_rhythmic[1, c("L5",  "L6",  "L7")]  <- 4      # T4
data_rhythmic[1, c("L9",  "L10", "L11")] <- 8      # T8
data_rhythmic[1, c("L13", "L14", "L15")] <- 12     # T12
data_rhythmic[1, c("L17", "L18", "L19")] <- 16     # T16
data_rhythmic[1, c("L21", "L22", "L23")] <- 20     # T20
data_rhythmic[1, c("L25", "L26", "L27")] <- 24     # T24
data_rhythmic[1, c("L29", "L30", "L31")] <- 28     # T28
data_rhythmic[1, c("L33", "L34", "L35")] <- 32     # T32
data_rhythmic[1, c("L37", "L38", "L39")] <- 36     # T36
data_rhythmic[1, c("L41", "L42", "L43")] <- 40     # T40
data_rhythmic[1, c("L45", "L46", "L47")] <- 44     # T44
data_rhythmic[1, c("L49", "L50", "L51")] <- 48     # T48

data_rhythmic <- data_rhythmic %>% mutate_if(is.character, as.numeric)

############################################################
# II. SEPARATE DATA BY REPLICATE
############################################################

rep1 <- data_rhythmic[, c("L1","L5","L9","L13","L17","L21","L25","L29","L33","L37","L41","L45","L49")]
rep2 <- data_rhythmic[, c("L2","L6","L10","L14","L18","L22","L26","L30","L34","L38","L42","L46","L50")]
rep3 <- data_rhythmic[, c("L3","L7","L11","L15","L19","L23","L27","L31","L35","L39","L43","L47","L51")]

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

rhythmic_TF  <- filter(TF, ID %in% cycling_TOT_leaf$AGI)
input.genes  <- unique(rhythmic_TF$ID)

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
            "outputs/data/link_list_dynGENIE3.txt",
            row.names = FALSE,
            quote     = FALSE,
            sep       = "\t")