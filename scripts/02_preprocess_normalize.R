#!/usr/bin/env Rscript
################################################################################
# STEP 02: PREPROCESS & NORMALIZE
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(SummarizedExperiment)
})

cat("\n=== PREPROCESSING ===\n\n")

mirna_raw <- readRDS("data/raw/mirna_raw.rds")
mrna_raw <- readRDS("data/raw/mrna_raw.rds")
clinical_raw <- readRDS("data/raw/clinical_raw.rds")

## Process miRNA
cat("Processing miRNA...\n")

if(is(mirna_raw, "SummarizedExperiment")) {
  mirna_data <- assay(mirna_raw, 1)
} else {
  mirna_data <- as.matrix(mirna_raw)
}

# Strip assay prefix (read_count_TCGA-... -> TCGA-...)
colnames(mirna_data) <- str_replace(colnames(mirna_data), "^.*_(TCGA-.+)$", "\\1")

# Remove duplicates (each sample appears 3x for 3 assay types)
mirna_data <- mirna_data[, !duplicated(colnames(mirna_data))]

# Filter primary tumors
primary_idx <- str_detect(colnames(mirna_data), "-01A-")
mirna_data <- mirna_data[, primary_idx]

# Handle pre-normalized data
if(any(mirna_data < 0, na.rm=TRUE) || max(mirna_data, na.rm=TRUE) < 1000) {
  keep <- apply(mirna_data, 1, var, na.rm=TRUE) > 0.01
  mirna_log2cpm <- mirna_data[keep, ]
} else {
  keep <- rowSums(cpm(mirna_data) > 1) >= ceiling(ncol(mirna_data)*0.1)
  mirna_log2cpm <- log2(cpm(mirna_data[keep,]) + 1)
}

cat(sprintf("  miRNA: %d x %d\n", nrow(mirna_log2cpm), ncol(mirna_log2cpm)))

## Process mRNA
cat("Processing mRNA...\n")

mrna_counts <- assay(mrna_raw, "unstranded")
primary_idx <- str_detect(colnames(mrna_counts), "-01A-")
mrna_counts <- mrna_counts[, primary_idx]

keep <- rowSums(cpm(mrna_counts) > 1) >= ceiling(ncol(mrna_counts)*0.1)
mrna_log2cpm <- log2(cpm(mrna_counts[keep,]) + 1)

cat(sprintf("  mRNA: %d x %d\n", nrow(mrna_log2cpm), ncol(mrna_log2cpm)))

## Match samples
cat("Matching samples...\n")
extract_bc <- function(x) str_sub(as.character(x), 1, 12)

mirna_bc <- extract_bc(colnames(mirna_log2cpm))
mrna_bc <- extract_bc(colnames(mrna_log2cpm))
common_bc <- intersect(mirna_bc, mrna_bc)

cat(sprintf("  Common: %d\n", length(common_bc)))

## Create balanced progression labels
cat("Creating progression labels...\n")

clinical_matched <- clinical_raw %>%
  mutate(patient_barcode = extract_bc(submitter_id)) %>%
  filter(patient_barcode %in% common_bc) %>%
  mutate(
    stage_clean = str_extract(ajcc_pathologic_stage, "Stage [IVX]+[AB]?"),
    progression_group = case_when(
      str_detect(stage_clean, "Stage I[^IV]|Stage IA|Stage IB|Stage IIA") ~ "early",
      str_detect(stage_clean, "Stage IIB|Stage III|Stage IV") ~ "late",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(progression_group)) %>%
  select(patient_barcode, barcode=submitter_id, stage_raw=ajcc_pathologic_stage, progression_group)

cat(sprintf("  Labeled: %d (Early: %d [%.1f%%], Late: %d [%.1f%%])\n\n",
            nrow(clinical_matched),
            sum(clinical_matched$progression_group == "early"),
            100*mean(clinical_matched$progression_group == "early"),
            sum(clinical_matched$progression_group == "late"),
            100*mean(clinical_matched$progression_group == "late")))

## Save
dir.create("data/processed", showWarnings=FALSE, recursive=TRUE)
saveRDS(mirna_log2cpm, "data/processed/mirna_log2cpm.rds")
saveRDS(mrna_log2cpm, "data/processed/mrna_log2cpm.rds")
saveRDS(clinical_matched, "data/processed/clinical_progression.rds")

cat("✓ Complete\n\n")
