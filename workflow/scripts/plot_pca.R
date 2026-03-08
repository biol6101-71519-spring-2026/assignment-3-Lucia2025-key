#!/usr/bin/env Rscript
#
# Generate PCA plot from count matrix
#
# This script is called by Snakemake and has access to:
# - snakemake@input: list of input files
# - snakemake@output: list of output files  
# - snakemake@params: parameters from the rule
# - snakemake@log: log file path

library(ggplot2)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

cat("Reading count matrix...\n")
counts <- read.table(snakemake@input$counts, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cat("Count matrix dimensions:", dim(counts), "\n")

cat("Reading sample metadata...\n")
metadata <- read.table(snakemake@input$metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Sample metadata:\n")
print(metadata)

cat("\nFiltering low-count genes...\n")
keep <- rowSums(counts >= 10) >= 3
counts_filtered <- counts[keep, , drop = FALSE]
cat("Genes after filtering:", nrow(counts_filtered), "\n")

cat("\nSelecting top variable genes...\n")
n_top <- min(500, nrow(counts_filtered))
gene_vars <- apply(counts_filtered, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_top]
counts_top <- counts_filtered[top_genes, , drop = FALSE]
cat("Using top", nrow(counts_top), "most variable genes\n")

cat("\nLog-transforming counts...\n")
counts_log <- log2(counts_top + 1)

cat("Removing zero-variance genes...\n")
gene_vars_log <- apply(counts_log, 1, var)
counts_log <- counts_log[gene_vars_log > 0, , drop = FALSE]
cat("Genes after removing zero-variance genes:", nrow(counts_log), "\n")

if (nrow(counts_log) < 2) {
  cat("No variable genes detected. Creating placeholder PCA plot.\n")
  pca_df <- data.frame(
    PC1 = rep(0, nrow(metadata)),
    PC2 = rep(0, nrow(metadata)),
    sample = metadata$sample,
    condition = metadata$condition
  )
  var_explained <- c(0, 0)
} else {
  cat("\nPerforming PCA...\n")
  pca_result <- prcomp(t(counts_log), scale. = FALSE, center = TRUE)
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

  pca_df <- as.data.frame(pca_result$x[, 1:2, drop = FALSE])
  pca_df$sample <- rownames(pca_df)
  pca_df$condition <- metadata$condition[match(pca_df$sample, metadata$sample)]
}

cat("\nGenerating PCA plot...\n")
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -1, size = 3) +
  labs(
    title = "PCA of RNA-Seq Samples",
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance)"),
    color = "Condition"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(snakemake@output[[1]], plot = p, width = 8, height = 6)
cat("PCA plot saved to", snakemake@output[[1]], "\n")

sink(type = "message")
sink(type = "output")
close(log_file)