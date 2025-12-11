#!/usr/bin/env Rscript
################################################################################
# STEP 07: EDA & DIFFERENTIAL EXPRESSION
# Robust version: ensure numeric, remove NAs, avoid integer overflow
################################################################################

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("DESeq2")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
})

cat("\n=== EDA & DE ANALYSIS ===\n\n")

mirna_raw <- readRDS("data/processed/mirna_log2cpm.rds")
clinical   <- readRDS("data/processed/clinical_integrated/clinical_ml_features.rds")

cat("Coercing miRNA matrix to numeric...\n")

# Coerce to numeric matrix safely
mirna <- as.matrix(mirna_raw)
mode(mirna) <- "numeric"

# Restore / set rownames
if (!is.null(rownames(mirna_raw))) {
  rownames(mirna) <- rownames(mirna_raw)
} else if (is.null(rownames(mirna))) {
  rownames(mirna) <- paste0("miRNA_", seq_len(nrow(mirna)))
}

cat(sprintf("  miRNA: %d features x %d samples\n", nrow(mirna), ncol(mirna)))

cat("  Range(log2CPM) before cleanup: ",
    paste(round(range(mirna, na.rm = TRUE), 3), collapse = " - "),
    "\n", sep = "")

# Replace non-finite values with NA, then drop problematic rows/columns
mirna[!is.finite(mirna)] <- NA

# Drop features with any NA
drop_rows <- rowSums(is.na(mirna)) > 0
if (any(drop_rows)) {
  cat(sprintf("  Dropping %d features with NA values\n", sum(drop_rows)))
  mirna <- mirna[!drop_rows, , drop = FALSE]
}

# Drop samples with any NA (should be rare)
drop_cols <- colSums(is.na(mirna)) > 0
if (any(drop_cols)) {
  cat(sprintf("  Dropping %d samples with NA values\n", sum(drop_cols)))
  mirna <- mirna[, !drop_cols, drop = FALSE]
}

cat("Computing PCA...\n")
pca <- prcomp(t(mirna), scale. = TRUE)
cat(sprintf("  PC1 explains %.1f%% variance\n", summary(pca)$importance[2, 1] * 100))
cat(sprintf("  PC2 explains %.1f%% variance\n\n", summary(pca)$importance[2, 2] * 100))

# Differential expression
cat("Preparing count matrix for DESeq2...\n")

# Back-transform log2CPM to pseudo-counts, avoid negatives / overflow
mirna_lin <- 2^mirna - 1
mirna_lin[mirna_lin < 0 | !is.finite(mirna_lin)] <- 0

int_max <- .Machine$integer.max
too_big <- mirna_lin > int_max
if (any(too_big)) {
  n_big <- sum(too_big)
  cat(sprintf("  WARNING: %d values > %d, capping to max integer\n",
              n_big, int_max))
  mirna_lin[too_big] <- int_max
}

mirna_counts <- round(mirna_lin)

if (anyNA(mirna_counts)) {
  stop("NA values remain in mirna_counts after cleanup. Check input matrix.")
}

# Match samples to clinical
extract_bc <- function(x) stringr::str_sub(as.character(x), 1, 12)

sample_bc <- extract_bc(colnames(mirna_counts))
clinical_matched <- clinical[match(sample_bc, clinical$patient_barcode), ]

keep_idx <- !is.na(clinical_matched$progression_group)
mirna_counts     <- mirna_counts[, keep_idx, drop = FALSE]
clinical_matched <- clinical_matched[keep_idx, , drop = FALSE]

cat(sprintf("  Matched %d samples with clinical data\n", sum(keep_idx)))

# DESeq2 dataset
coldata <- data.frame(
  condition = factor(clinical_matched$progression_group,
                     levels = c("early", "late")),
  row.names = colnames(mirna_counts)
)

dds <- DESeqDataSetFromMatrix(
  countData = mirna_counts,
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds, quiet = TRUE)

results_de <- as.data.frame(results(dds, contrast = c("condition", "late", "early"))) %>%
  tibble::rownames_to_column("feature") %>%
  arrange(padj)

sig <- results_de %>%
  filter(!is.na(padj),
         padj < 0.05,
         abs(log2FoldChange) > 1)

cat(sprintf("  Significant: %d features (padj<0.05, |LFC|>1)\n\n", nrow(sig)))

# Save outputs
dir.create("results/eda", showWarnings = FALSE, recursive = TRUE)

readr::write_csv(results_de, "results/eda/de_results.csv")
readr::write_csv(sig,        "results/eda/de_significant.csv")
saveRDS(list(pca = pca, de = results_de), "results/eda/eda_summary.rds")

cat("=== EDA & DE ANALYSIS COMPLETE ===\n\n")

if (nrow(sig) > 0) {
  cat("Top 10 DE features:\n")
  print(head(sig %>% dplyr::select(feature, log2FoldChange, padj), 10))
  cat("\n")
}

cat("Files saved:\n")
cat("  - results/eda/de_results.csv\n")
cat("  - results/eda/de_significant.csv\n")
cat("  - results/eda/eda_summary.rds\n\n")

cat("Next: Rscript scripts/08_ml_benchmarking.R\n\n")
