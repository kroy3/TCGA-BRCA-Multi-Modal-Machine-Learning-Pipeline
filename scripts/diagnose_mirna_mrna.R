#!/usr/bin/env Rscript
################################################################################
# DIAGNOSTIC: Check miRNA/mRNA data structure
################################################################################

suppressPackageStartupMessages(library(tidyverse))

cat("\n=== DIAGNOSTIC: miRNA-mRNA MATCHING ===\n\n")

mirna <- readRDS("data/processed/mirna_log2cpm.rds")
mrna <- readRDS("data/processed/mrna_log2cpm.rds")

cat("1. DATA DIMENSIONS:\n")
cat(sprintf("   miRNA: %d features x %d samples\n", nrow(mirna), ncol(mirna)))
cat(sprintf("   mRNA:  %d features x %d samples\n\n", nrow(mrna), ncol(mrna)))

cat("2. SAMPLE NAMES (first 5):\n")
cat("   miRNA:\n")
print(colnames(mirna)[1:5])
cat("\n   mRNA:\n")
print(colnames(mrna)[1:5])

cat("\n3. BARCODE EXTRACTION TESTS:\n")

# Test different extraction lengths
for(len in c(12, 15, 16, 19)) {
  extract_bc <- function(x) str_sub(as.character(x), 1, len)
  
  mirna_bc <- extract_bc(colnames(mirna))
  mrna_bc <- extract_bc(colnames(mrna))
  common_bc <- intersect(mirna_bc, mrna_bc)
  
  cat(sprintf("   %d-char: %d common samples\n", len, length(common_bc)))
}

cat("\n4. DIRECT COLUMN NAME MATCH:\n")
direct_match <- intersect(colnames(mirna), colnames(mrna))
cat(sprintf("   Direct intersect: %d samples\n", length(direct_match)))

cat("\n5. CHECK DATA STRUCTURE:\n")
cat(sprintf("   miRNA class: %s\n", class(mirna)[1]))
cat(sprintf("   mRNA class:  %s\n", class(mrna)[1]))
cat(sprintf("   miRNA is matrix: %s\n", is.matrix(mirna)))
cat(sprintf("   mRNA is matrix:  %s\n", is.matrix(mrna)))

cat("\n6. CHECK VARIANCE COMPUTATION:\n")
tryCatch({
  mirna_var <- apply(mirna, 1, var, na.rm=TRUE)
  cat(sprintf("   miRNA variance: OK (%d values)\n", length(mirna_var)))
  cat(sprintf("   Range: %.4f to %.4f\n", min(mirna_var), max(mirna_var)))
  cat(sprintf("   Non-zero: %d\n", sum(mirna_var > 0)))
}, error = function(e) {
  cat(sprintf("   miRNA variance: ERROR - %s\n", e$message))
})

cat("\n7. BEST BARCODE LENGTH:\n")
extract_bc <- function(x) str_sub(as.character(x), 1, 16)
mirna_bc <- extract_bc(colnames(mirna))
mrna_bc <- extract_bc(colnames(mrna))
common_bc <- intersect(mirna_bc, mrna_bc)

cat(sprintf("   Using 16-char: %d common samples\n", length(common_bc)))
cat("   miRNA barcodes (first 3):\n")
print(mirna_bc[1:3])
cat("   mRNA barcodes (first 3):\n")
print(mrna_bc[1:3])

cat("\n8. TEST SUBSET:\n")
if(length(common_bc) > 0) {
  mirna_idx <- match(common_bc, mirna_bc)
  mrna_idx <- match(common_bc, mrna_bc)
  
  cat(sprintf("   mirna_idx length: %d\n", length(mirna_idx)))
  cat(sprintf("   mrna_idx length: %d\n", length(mrna_idx)))
  cat(sprintf("   NA in mirna_idx: %d\n", sum(is.na(mirna_idx))))
  cat(sprintf("   NA in mrna_idx: %d\n", sum(is.na(mrna_idx))))
  
  mirna_matched <- mirna[, mirna_idx]
  mrna_matched <- mrna[, mrna_idx]
  
  cat(sprintf("   Matched miRNA: %d x %d\n", nrow(mirna_matched), ncol(mirna_matched)))
  cat(sprintf("   Matched mRNA:  %d x %d\n\n", nrow(mrna_matched), ncol(mrna_matched)))
  
  # Test variance on matched data
  cat("9. VARIANCE ON MATCHED DATA:\n")
  mirna_var <- apply(mirna_matched, 1, var, na.rm=TRUE)
  cat(sprintf("   miRNA variance computed: %d values\n", length(mirna_var)))
  cat(sprintf("   Top variance: %.4f\n", max(mirna_var)))
  cat(sprintf("   Non-zero variance: %d\n", sum(mirna_var > 0)))
  
  top50_mirna <- names(sort(mirna_var, decreasing=TRUE)[1:min(50, sum(mirna_var > 0))])
  cat(sprintf("   Top variable miRNAs selected: %d\n", length(top50_mirna)))
  
  if(length(top50_mirna) > 0) {
    cat("\n   First top miRNA:\n")
    cat(sprintf("     Name: %s\n", top50_mirna[1]))
    mirna_expr <- as.numeric(mirna_matched[top50_mirna[1], ])
    cat(sprintf("     Expression length: %d\n", length(mirna_expr)))
    cat(sprintf("     Range: %.4f to %.4f\n", min(mirna_expr, na.rm=TRUE), max(mirna_expr, na.rm=TRUE)))
  }
}

cat("\n=== DIAGNOSTIC COMPLETE ===\n\n")
