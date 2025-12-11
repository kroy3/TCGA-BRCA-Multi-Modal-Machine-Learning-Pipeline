#!/usr/bin/env Rscript
################################################################################
# STEP 01: DOWNLOAD TCGA-BRCA DATA
################################################################################

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

cat("\n=== DOWNLOADING TCGA-BRCA DATA ===\n\n")

dir.create("data/raw", recursive=TRUE, showWarnings=FALSE)

## 1. mRNA - Keep only unstranded assay to save memory
cat("1. Downloading mRNA (30-45 min)...\n")
query_mrna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_mrna, method="api", files.per.chunk=50)
mrna_full <- GDCprepare(query_mrna)

# Extract only unstranded counts
mrna_data <- SummarizedExperiment(
  assays = list(unstranded = assay(mrna_full, "unstranded")),
  colData = colData(mrna_full),
  rowData = rowData(mrna_full)
)
rm(mrna_full); gc()

saveRDS(mrna_data, "data/raw/mrna_raw.rds")
cat(sprintf("   ✓ %d genes x %d samples\n\n", nrow(mrna_data), ncol(mrna_data)))

## 2. miRNA
cat("2. Downloading miRNA (15-30 min)...\n")
query_mirna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification"
)

GDCdownload(query_mirna, method="api", files.per.chunk=50, directory="GDCdata_miRNA")
mirna_data <- GDCprepare(query_mirna, directory="GDCdata_miRNA")
saveRDS(mirna_data, "data/raw/mirna_raw.rds")
cat(sprintf("   ✓ %d features x %d samples\n\n", nrow(mirna_data), ncol(mirna_data)))

## 3. Clinical
cat("3. Downloading clinical data...\n")
clinical_data <- GDCquery_clinic("TCGA-BRCA", "clinical")
saveRDS(clinical_data, "data/raw/clinical_raw.rds")
cat(sprintf("   ✓ %d patients\n\n", nrow(clinical_data)))

cat("=== DOWNLOAD COMPLETE ===\n")
cat("Next: Rscript scripts/02_preprocess_normalize.R\n\n")
