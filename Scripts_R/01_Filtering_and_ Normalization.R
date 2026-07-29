

############################################################
# Script: 01_filtering_normalization_DESeq2.R
# Project: Circadian_gene_expression_ analysis in Pisum Sativum  
#
# Description:
# Filtering of low-expressed genes, normalization of RNA-seq
# count data using DESeq2 variance stabilizing transformation (VST),
# and generation of normalized expression matrices for downstream
# circadian and network analyses.
#
# Input:
# - Raw gene count matrix
# - Sample metadata file
#
# Output:
# - Normalized expression matrix
# - Leaf-specific expression matrix
# - Seed-specific expression matrix
#
# Author: Fatoumata KAMISSOKO 
############################################################


## Load libraries -------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(ade4)
library(factoextra)


## Define paths ---------------------------------------------------------------

count_file <- "data/raw/salmon.merged.gene_counts.tsv"
metadata_file <- "data/raw/Metadata.txt"

output_dir <- "outputs/data"


## Import count matrix and metadata --------------------------------------------

data <- read.table(
  count_file,
  header = TRUE,
  sep = "\t"
)

metadata <- read.table(
  metadata_file,
  header = TRUE,
  sep = "\t"
)


## Clean column names ----------------------------------------------------------

names(data) <- gsub(
  pattern = "read_count_",
  replacement = "",
  x = names(data)
)


## Set row names --------------------------------------------------------------

row.names(metadata) <- metadata$SampleName
row.names(data) <- data$gene_id


## Convert metadata variables -----------------------------------------------

metadata$Time <- as.factor(metadata$Time)
metadata$Organ <- as.factor(metadata$Organ)


## Check sample consistency ----------------------------------------------------

all(rownames(metadata) %in% colnames(data))

# Reorder count matrix according to metadata
data <- data[, metadata$SampleName, drop = FALSE]

all(rownames(metadata) == colnames(data))


## Convert fractional counts to integers --------------------------------------

data.round.counts <- data %>%
  mutate_if(is.numeric, as.integer)


## Filter low expressed genes --------------------------------------------------

# Keep genes with >=10 reads in at least 50% of samples

keep <- rowSums(data.round.counts >= 10) >= 0.5*ncol(data.round.counts)

data.filtered <- data.round.counts[keep, ]

cat(
  "Number of retained genes:",
  sum(keep),
  "\n"
)


## DESeq2 normalization using VST ---------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = data.filtered,
  colData = metadata,
  design = ~ Time
)


vst_data <- vst(dds)

Expressed.vst <- assay(vst_data)

Expressed.vst <- data.frame(Expression = Expressed.vst)

Expressed.vst$ID <- row.names(Expressed.vst)

Expressed.vst <- Expressed.vst[, c("ID", row.names(metadata))]


## Export normalized expression matrix -----------------------------------------

write.table(
  Expressed.vst,
  file.path(output_dir, "Expressed.vst.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## Separate leaves and seeds ---------------------------------------------------

samples_seeds <- metadata$SampleName[
  metadata$Organ == "S"
]

samples_leaves <- metadata$SampleName[
  metadata$Organ == "L"
]


expression.seeds <- Expressed.vst[, c("ID", samples_seeds)]

expression.leaves <- Expressed.vst[, c("ID", samples_leaves)]


write.table(
  expression.seeds,
  file.path(output_dir, "expression.seeds.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  expression.leaves,
  file.path(output_dir, "expression.leaves.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

