#!/usr/bin/env Rscript
#
# Prepare DESeq2 dataset for differential expression analysis
#
# This script creates a DESeq2 dataset object and performs normalization.
#
# This script is called by Snakemake and has access to:
# - snakemake@input: list of input files
# - snakemake@output: list of output files
# - snakemake@params: parameters from the rule
# - snakemake@log: log file path

# Redirect output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
})

# Read count matrix
cat("Reading count matrix...\n")
counts <- read.table(snakemake@input$counts, header = TRUE, row.names = 1, sep = "\t")
cat("Count matrix dimensions:", dim(counts), "\n")
cat("Total reads per sample:\n")
print(colSums(counts, na.rm = TRUE))

# Read sample metadata
cat("\nReading sample metadata...\n")
metadata <- read.table(snakemake@input$metadata, header = TRUE, sep = "\t", row.names = 1)
cat("Sample metadata:\n")
print(metadata)

# Ensure sample order matches
counts <- counts[, rownames(metadata), drop = FALSE]
cat("\nSample order verified.\n")

# Replace NA values with 0
cat("\nReplacing NA values with 0...\n")
counts[is.na(counts)] <- 0
cat("Number of NA values after cleanup:", sum(is.na(counts)), "\n")

# Force counts to be a numeric matrix with preserved row/column names
counts_mat <- as.matrix(counts)
mode(counts_mat) <- "numeric"
storage.mode(counts_mat) <- "integer"
rownames(counts_mat) <- rownames(counts)
colnames(counts_mat) <- colnames(counts)

cat("Count matrix dimensions after cleanup:", dim(counts_mat), "\n")

# Filter low-count genes
cat("\nFiltering low-count genes...\n")
min_count <- if (!is.null(snakemake@params$min_count)) snakemake@params$min_count else 10
cat("Minimum count threshold:", min_count, "\n")

keep <- rowSums(counts_mat >= min_count, na.rm = TRUE) >= 3
counts_filtered <- counts_mat[keep, , drop = FALSE]

cat("Genes before filtering:", nrow(counts_mat), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")

# Create DESeq2 dataset
cat("\nCreating DESeq2 dataset...\n")
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ condition
)

cat("DESeq2 dataset created successfully.\n")
cat("Design formula: ~ condition\n")
cat("Samples:", ncol(dds), "\n")
cat("Genes:", nrow(dds), "\n")

# Run DESeq2 normalization
cat("\nRunning DESeq2 normalization...\n")
dds <- estimateSizeFactors(dds)

cat("Size factors:\n")
print(sizeFactors(dds))

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

cat("\nNormalized counts summary:\n")
cat("Total normalized reads per sample:\n")
print(colSums(normalized_counts))

# Save DESeq2 dataset object
cat("\nSaving DESeq2 dataset...\n")
saveRDS(dds, file = snakemake@output$dds)
cat("DESeq2 dataset saved to", snakemake@output$dds, "\n")

# Save normalized counts
write.table(
  normalized_counts,
  file = snakemake@output$normalized,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)
cat("Normalized counts saved to", snakemake@output$normalized, "\n")

cat("\n=== DESeq2 dataset ready for differential expression analysis ===\n")

# Close log connections
sink(type = "message")
sink(type = "output")
close(log_file)