#!/usr/bin/env Rscript
################################################################################
# STEP 05: CLINICAL FEATURE INTEGRATION
################################################################################

if(!require(TCGAbiolinks, quietly=TRUE)) {
  if(!require(BiocManager, quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install("TCGAbiolinks")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(TCGAbiolinks)
})

cat("\n=== CLINICAL INTEGRATION ===\n\n")

clinical <- readRDS("data/processed/clinical_progression.rds")

cat("Downloading TCGA clinical data...\n")
tcga_clin <- GDCquery_clinic("TCGA-BRCA", "clinical")

extract_bc <- function(x) str_sub(as.character(x), 1, 12)
clinical$patient_barcode <- extract_bc(clinical$barcode)
tcga_clin$patient_barcode <- extract_bc(tcga_clin$submitter_id)

# Merge and encode
clinical_ml <- clinical %>%
  left_join(tcga_clin %>% select(patient_barcode, age_at_index, 
            ajcc_pathologic_t, ajcc_pathologic_n, ajcc_pathologic_m,
            race, gender, vital_status), by="patient_barcode") %>%
  mutate(
    age = as.numeric(age_at_index),
    t_stage = as.numeric(str_extract(ajcc_pathologic_t, "[0-4]")),
    n_stage = as.numeric(str_extract(ajcc_pathologic_n, "[0-3]")),
    has_metastasis = as.numeric(str_detect(ajcc_pathologic_m, "M1")),
    is_female = as.numeric(str_detect(tolower(gender), "female")),
    race_white = as.numeric(str_detect(tolower(race), "white")),
    died = as.numeric(str_detect(tolower(vital_status), "dead")),
    # Impute missing
    age = ifelse(is.na(age), median(age, na.rm=TRUE), age),
    t_stage = ifelse(is.na(t_stage), 2, t_stage),
    n_stage = ifelse(is.na(n_stage), 0, n_stage)
  ) %>%
  select(patient_barcode, age, t_stage, n_stage, has_metastasis, 
         is_female, race_white, died, progression_group)

dir.create("data/processed/clinical_integrated", showWarnings=FALSE, recursive=TRUE)
saveRDS(clinical_ml, "data/processed/clinical_integrated/clinical_ml_features.rds")

cat(sprintf("✓ %d samples with 7 clinical features\n\n", nrow(clinical_ml)))
