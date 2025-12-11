#!/usr/bin/env Rscript
################################################################################
# DEEP DIAGNOSTIC: Check RDS file structure
################################################################################

suppressPackageStartupMessages(library(tidyverse))

cat("\n=== DEEP RDS DIAGNOSTIC ===\n\n")

cat("1. LOADING FILES (RAW):\n")
mirna_raw <- readRDS("data/processed/mirna_log2cpm.rds")
mrna_raw <- readRDS("data/processed/mrna_log2cpm.rds")

cat(sprintf("   mirna_raw class: %s\n", paste(class(mirna_raw), collapse=", ")))
cat(sprintf("   mrna_raw class: %s\n", paste(class(mrna_raw), collapse=", ")))

cat("\n2. STRUCTURE:\n")
cat("   miRNA structure:\n")
str(mirna_raw)

cat("\n   mRNA structure:\n")
str(mrna_raw)

cat("\n3. DIMENSIONS:\n")
cat(sprintf("   miRNA: %s\n", paste(dim(mirna_raw), collapse=" x ")))
cat(sprintf("   mRNA: %s\n", paste(dim(mrna_raw), collapse=" x ")))

cat("\n4. ROWNAMES:\n")
mirna_rn <- rownames(mirna_raw)
mrna_rn <- rownames(mrna_raw)

cat(sprintf("   miRNA rownames: %s\n", 
            ifelse(is.null(mirna_rn), "NULL", paste0(length(mirna_rn), " names"))))
cat(sprintf("   mRNA rownames: %s\n", 
            ifelse(is.null(mrna_rn), "NULL", paste0(length(mrna_rn), " names"))))

if(!is.null(mirna_rn) && length(mirna_rn) > 0) {
  cat("   miRNA rownames (first 5):\n")
  print(head(mirna_rn, 5))
}

if(!is.null(mrna_rn) && length(mrna_rn) > 0) {
  cat("   mRNA rownames (first 5):\n")
  print(head(mrna_rn, 5))
}

cat("\n5. COLNAMES:\n")
mirna_cn <- colnames(mirna_raw)
mrna_cn <- colnames(mrna_raw)

cat(sprintf("   miRNA colnames: %s\n", 
            ifelse(is.null(mirna_cn), "NULL", paste0(length(mirna_cn), " names"))))
cat(sprintf("   mRNA colnames: %s\n", 
            ifelse(is.null(mrna_cn), "NULL", paste0(length(mrna_cn), " names"))))

if(!is.null(mirna_cn) && length(mirna_cn) > 0) {
  cat("   miRNA colnames (first 5):\n")
  print(head(mirna_cn, 5))
}

if(!is.null(mrna_cn) && length(mrna_cn) > 0) {
  cat("   mRNA colnames (first 5):\n")
  print(head(mrna_cn, 5))
}

cat("\n6. DIMNAMES:\n")
cat("   miRNA dimnames:\n")
print(str(dimnames(mirna_raw)))

cat("\n   mRNA dimnames:\n")
print(str(dimnames(mrna_raw)))

cat("\n7. ATTRIBUTES:\n")
cat("   miRNA attributes:\n")
print(names(attributes(mirna_raw)))

cat("\n   mRNA attributes:\n")
print(names(attributes(mrna_raw)))

cat("\n8. IF DATA.FRAME:\n")
if(is.data.frame(mirna_raw)) {
  cat("   miRNA is data.frame\n")
  cat("   Column names:\n")
  print(head(names(mirna_raw), 10))
}

if(is.data.frame(mrna_raw)) {
  cat("   mRNA is data.frame\n")
  cat("   Column names:\n")
  print(head(names(mrna_raw), 10))
}

cat("\n9. ACCESSING DATA:\n")
cat("   Trying mirna_raw[1:3, 1:3]:\n")
print(mirna_raw[1:3, 1:3])

cat("\n10. CHECKING PREPROCESSING SCRIPT:\n")
if(file.exists("data/processed/clinical_progression.rds")) {
  clinical <- readRDS("data/processed/clinical_progression.rds")
  cat(sprintf("   Clinical samples: %d\n", nrow(clinical)))
  cat("   Clinical structure:\n")
  print(str(clinical))
}

cat("\n=== END DIAGNOSTIC ===\n\n")
