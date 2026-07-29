
############################################################
# Script: 02_MetaCycle_rhythmic_detection.R
# Project: Circadian gene expression analysis 
#
# Description:
# Identification of rhythmic genes from normalized transcriptomic
# expression matrices using MetaCycle. Rhythmic gene detection
# is performed independently in leaves and seeds using multiple
# algorithms integrated in MetaCycle.
#
# Input:
# - Normalized expression matrix for leaves
# - Normalized expression matrix for seeds
#
# Output:
# - Complete MetaCycle results for leaves and seeds
# - Lists of significantly rhythmic genes
#
# Author: Fatoumata KAMISSOKO 
############################################################


## Load libraries -------------------------------------------------------------

library(MetaCycle)


## Define paths ---------------------------------------------------------------

leaf_expression_file <- "outputs/data/expression.leaves.txt"
seed_expression_file <- "outputs/data/expression.seeds.txt"

output_dir <- "outputs/data"


## Load normalized expression data ---------------------------------------------

expression.leaves <- read.table(
  leaf_expression_file,
  header = TRUE,
  sep = "\t"
)


expression.seeds <- read.table(
  seed_expression_file,
  header = TRUE,
  sep = "\t"
)



## Define experimental time points ---------------------------------------------

# 48 hours experiment
# Sampling every 4 hours
# 3 biological replicates per time point

time_points <- rep(seq(0, 48, by = 4), each = 3)



## MetaCycle analysis in leaves -----------------------------------------------

Leaf_MetaCycle <- meta2d(
  infile = "txt",
  filestyle = "txt",
  minper = 20,
  maxper = 28,
  timepoints = time_points,
  outputFile = FALSE,
  adjustPhase = "notAdjusted",
  combinePvalue = "fisher",
  inDF = expression.leaves
)


Leaf_MetaCycle_results <- Leaf_MetaCycle$meta


# Significant rhythmic genes (BH adjusted p-value < 0.05)

Cycling_leaf <- Leaf_MetaCycle_results[
  Leaf_MetaCycle_results$meta2d_BH.Q < 0.05,
]



## MetaCycle analysis in seeds ------------------------------------------------

Seed_MetaCycle <- meta2d(
  infile = "txt",
  filestyle = "txt",
  minper = 20,
  maxper = 28,
  timepoints = time_points,
  outputFile = FALSE,
  adjustPhase = "notAdjusted",
  combinePValue = "fisher",
  inDF = expression.seeds
)


Seed_MetaCycle_results <- Seed_MetaCycle$meta


# Significant rhythmic genes

Cycling_seed <- Seed_MetaCycle_results[
  Seed_MetaCycle_results$meta2d_BH.Q < 0.05,
]



## Export results -------------------------------------------------------------

write.table(
  Leaf_MetaCycle_results,
  file.path(output_dir, "Leaf_MetaCycle_results.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  Seed_MetaCycle_results,
  file.path(output_dir, "Seed_MetaCycle_results.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



write.table(
  Cycling_leaf,
  file.path(output_dir, "Cycling_leaf.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  Cycling_seed,
  file.path(output_dir, "Cycling_seeds.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



## Summary --------------------------------------------------------------------

cat(
  "Number of rhythmic genes in leaves:",
  nrow(Cycling_leaf),
  "\n"
)

cat(
  "Number of rhythmic genes in seeds:",
  nrow(Cycling_seed),
  "\n"
)
