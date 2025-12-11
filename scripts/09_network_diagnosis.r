#!/usr/bin/env Rscript
################################################################################
# DIAGNOSTIC SCRIPT for 09_network_analysis
# Run with: Rscript scripts/09_network_diagnosis.R
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("\n=== NETWORK DIAGNOSTIC ===\n\n")

mirna_raw <- readRDS("data/processed/mirna_log2cpm.rds")
mrna_raw  <- readRDS("data/processed/mrna_log2cpm.rds")

cat("--- Basic structure ---\n")
cat("miRNA class: ", paste(class(mirna_raw), collapse = ","), "\n")
cat("mRNA  class: ", paste(class(mrna_raw),  collapse = ","), "\n")
cat(sprintf("miRNA dim: %d x %d\n", nrow(mirna_raw), ncol(mirna_raw)))
cat(sprintf("mRNA  dim: %d x %d\n\n", nrow(mrna_raw), ncol(mrna_raw)))

cat("First 5 miRNA column names:\n")
print(head(colnames(mirna_raw), 5))
cat("\nFirst 5 mRNA column names:\n")
print(head(colnames(mrna_raw), 5))
cat("\n")

# Coerce to numeric matrices ---------------------------------------------------
mirna <- as.matrix(mirna_raw)
mrna  <- as.matrix(mrna_raw)
mode(mirna) <- "numeric"
mode(mrna)  <- "numeric"

cat("--- Numeric checks ---\n")
cat("is.numeric(mirna[1,1]): ", is.numeric(mirna[1, 1]), "\n")
cat("is.numeric(mrna[1,1]) : ", is.numeric(mrna[1, 1]),  "\n")
cat("Any NA in miRNA? ", anyNA(mirna), "\n")
cat("Any NA in mRNA?  ", anyNA(mrna),  "\n\n")

cat("--- Rowname status ---\n")
cat("miRNA has rownames? ", !is.null(rownames(mirna)), "\n")
if (!is.null(rownames(mirna))) {
  cat("First 5 miRNA rownames:\n")
  print(head(rownames(mirna), 5))
}
cat("mRNA has rownames?  ", !is.null(rownames(mrna)), "\n")
if (!is.null(rownames(mrna))) {
  cat("First 5 mRNA rownames:\n")
  print(head(rownames(mrna), 5))
}
cat("\n")

# Align samples as in network script ------------------------------------------
cat("--- Sample alignment ---\n")
common <- intersect(colnames(mirna), colnames(mrna))
cat("Exact overlap length: ", length(common), "\n")

if (length(common) == 0) {
  cat("Trying 12-char TCGA barcodes...\n")
  extract_bc <- function(x) stringr::str_sub(as.character(x), 1, 12)
  mirna_bc <- extract_bc(colnames(mirna))
  mrna_bc  <- extract_bc(colnames(mrna))
  common_bc <- intersect(mirna_bc, mrna_bc)
  cat("Barcode overlap length: ", length(common_bc), "\n")
  if (length(common_bc) > 0) {
    common_bc <- sort(common_bc)
    mirna_idx <- match(common_bc, mirna_bc)
    mrna_idx  <- match(common_bc, mrna_bc)
    mirna <- mirna[, mirna_idx, drop = FALSE]
    mrna  <- mrna[, mrna_idx,  drop = FALSE]
    colnames(mirna) <- common_bc
    colnames(mrna)  <- common_bc
  }
} else {
  mirna <- mirna[, common, drop = FALSE]
  mrna  <- mrna[, common, drop = FALSE]
}

cat(sprintf("Aligned dims - miRNA: %d x %d\n", nrow(mirna), ncol(mirna)))
cat(sprintf("Aligned dims - mRNA : %d x %d\n\n", nrow(mrna),  ncol(mrna)))

# Variance diagnostics ---------------------------------------------------------
cat("--- Variance diagnostics ---\n")
var_mirna <- apply(mirna, 1, var, na.rm = TRUE)
var_mrna  <- apply(mrna, 1, var, na.rm = TRUE)

cat("miRNA var summary:\n")
print(summary(var_mirna))
cat("\n# miRNA with var>0: ", sum(is.finite(var_mirna) & var_mirna > 0), "\n")
cat("# miRNA with var==0 or NA: ", sum(!is.finite(var_mirna) | var_mirna <= 0), "\n\n")

cat("mRNA var summary:\n")
print(summary(var_mrna))
cat("\n# mRNA with var>0: ", sum(is.finite(var_mrna) & var_mrna > 0), "\n")
cat("# mRNA with var==0 or NA: ", sum(!is.finite(var_mrna) | var_mrna <= 0), "\n\n")

# Show how top features are being selected ------------------------------------
cat("--- Top feature selection test ---\n")

n_top_mirna <- min(20, nrow(mirna))
n_top_mrna  <- min(500, nrow(mrna))

cat("n_top_mirna: ", n_top_mirna, "\n")
cat("n_top_mrna : ", n_top_mrna,  "\n\n")

# 1) Using names(sort(...)) like the original script
vec_mirna_var <- sort(var_mirna, decreasing = TRUE)
cat("length(vec_mirna_var): ", length(vec_mirna_var), "\n")
cat("any(is.na(vec_mirna_var))? ", any(is.na(vec_mirna_var)), "\n")

candidate_by_names <- names(vec_mirna_var)[seq_len(min(5, length(vec_mirna_var)))]
cat("First 5 names(vec_mirna_var):\n")
print(candidate_by_names)

if (n_top_mirna > 0) {
  top_by_names <- names(vec_mirna_var)[seq_len(n_top_mirna)]
  cat("length(top_by_names): ", length(top_by_names), "\n")
}

# 2) Using rownames directly
if (!is.null(rownames(mirna))) {
  cat("First 5 rownames(mirna):\n")
  print(head(rownames(mirna), 5))
  top_by_rownames <- rownames(mirna)[seq_len(n_top_mirna)]
  cat("length(top_by_rownames): ", length(top_by_rownames), "\n")
}

cat("\n=== END DIAGNOSTIC ===\n")
