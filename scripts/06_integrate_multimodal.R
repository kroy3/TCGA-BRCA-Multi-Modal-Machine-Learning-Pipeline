#!/usr/bin/env Rscript
################################################################################
# STEP 06: MULTI-MODAL INTEGRATION
# Fixed: Use 12-char barcodes to match clinical data
################################################################################

if(!require(caret, quietly=TRUE)) install.packages("caret")

suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
})

cat("\n=== MULTI-MODAL INTEGRATION ===\n\n")

# Load all data
clinical <- readRDS("data/processed/clinical_integrated/clinical_ml_features.rds")
pathways_raw <- readRDS("data/processed/pathway_features/gsva_pathway_scores.rds")
mirna_raw <- readRDS("data/processed/mirna_log2cpm.rds")
interactions_raw <- readRDS("data/processed/mirna_target_features/mirna_target_interactions.rds")

cat("Loaded data:\n")
cat(sprintf("  Clinical: %d x %d\n", nrow(clinical), ncol(clinical)))
cat(sprintf("  Pathways: %d x %d\n", nrow(pathways_raw), ncol(pathways_raw)))
cat(sprintf("  miRNA: %d x %d\n", nrow(mirna_raw), ncol(mirna_raw)))
cat(sprintf("  Interactions: %d x %d\n\n", nrow(interactions_raw), ncol(interactions_raw)))

cat("Checking clinical barcodes:\n")
cat(sprintf("  First clinical barcode: %s\n", clinical$patient_barcode[1]))
cat(sprintf("  Barcode length: %d\n\n", nchar(clinical$patient_barcode[1])))

# Convert all to numeric matrices
cat("Converting to numeric...\n")

pathways <- apply(pathways_raw, 2, as.numeric)
rownames(pathways) <- rownames(pathways_raw)

mirna <- apply(mirna_raw, 2, as.numeric)
if(!is.null(rownames(mirna_raw))) {
  rownames(mirna) <- rownames(mirna_raw)
} else {
  rownames(mirna) <- paste0("miRNA_", 1:nrow(mirna))
}

interactions <- apply(interactions_raw, 2, as.numeric)
rownames(interactions) <- rownames(interactions_raw)

cat("  All data converted to numeric\n\n")

# CRITICAL: Extract 12-char barcodes to match clinical
extract_bc <- function(x) str_sub(as.character(x), 1, 12)

# View-specific selection and scaling
cat("Selecting and scaling features...\n")

select_top <- function(mat, n, name) {
  mat_num <- apply(mat, 2, as.numeric)
  rownames(mat_num) <- rownames(mat)
  
  vars <- apply(mat_num, 1, var, na.rm=TRUE)
  vars[is.na(vars)] <- 0
  
  n_actual <- min(n, sum(vars > 0))
  top_idx <- order(vars, decreasing=TRUE)[1:n_actual]
  mat_selected <- mat_num[top_idx, , drop=FALSE]
  
  mat_scaled <- t(scale(t(mat_selected)))
  rownames(mat_scaled) <- rownames(mat_selected)
  
  cat(sprintf("  %s: selected %d / %d features\n", name, nrow(mat_scaled), nrow(mat)))
  
  return(mat_scaled)
}

pathway_sel <- select_top(pathways, 150, "Pathways")
mirna_sel <- select_top(mirna, 100, "miRNA")
interact_sel <- t(scale(t(interactions)))
rownames(interact_sel) <- rownames(interactions)

cat("\n")

# Add prefixes
rownames(pathway_sel) <- paste0("PATHWAY_", rownames(pathway_sel))
rownames(mirna_sel) <- paste0("miRNA_", rownames(mirna_sel))
rownames(interact_sel) <- paste0("INTERACT_", rownames(interact_sel))

# Convert to data frames with 12-char barcodes
cat("Creating data frames with patient barcodes...\n")

pathway_df <- as.data.frame(t(pathway_sel))
pathway_df$patient_barcode <- extract_bc(rownames(pathway_df))

mirna_df <- as.data.frame(t(mirna_sel))
mirna_df$patient_barcode <- extract_bc(rownames(mirna_df))

interact_df <- as.data.frame(t(interact_sel))
interact_df$patient_barcode <- extract_bc(rownames(interact_df))

cat(sprintf("  Pathway barcodes: %d unique (first: %s)\n", 
            length(unique(pathway_df$patient_barcode)),
            pathway_df$patient_barcode[1]))
cat(sprintf("  miRNA barcodes: %d unique (first: %s)\n", 
            length(unique(mirna_df$patient_barcode)),
            mirna_df$patient_barcode[1]))
cat(sprintf("  Interaction barcodes: %d unique (first: %s)\n", 
            length(unique(interact_df$patient_barcode)),
            interact_df$patient_barcode[1]))

# Check overlap before merging
cat("\nChecking barcode overlap:\n")
cat(sprintf("  Clinical & Pathway: %d common\n", 
            length(intersect(clinical$patient_barcode, pathway_df$patient_barcode))))
cat(sprintf("  Clinical & miRNA: %d common\n", 
            length(intersect(clinical$patient_barcode, mirna_df$patient_barcode))))
cat(sprintf("  Clinical & Interactions: %d common\n\n", 
            length(intersect(clinical$patient_barcode, interact_df$patient_barcode))))

# Merge all
cat("Merging all data types...\n")

integrated <- clinical %>%
  inner_join(pathway_df, by="patient_barcode") %>%
  inner_join(mirna_df, by="patient_barcode") %>%
  inner_join(interact_df, by="patient_barcode") %>%
  filter(!is.na(progression_group))

cat(sprintf("Integrated: %d samples x %d features\n\n", 
            nrow(integrated), ncol(integrated)-2))

if(nrow(integrated) == 0) {
  cat("ERROR: No samples after integration\n")
  cat("Debugging information saved\n")
  quit(status=1)
}

# Create stratified splits
cat("Creating train/val/test splits...\n")
set.seed(42)

train_idx <- createDataPartition(integrated$progression_group, p=0.70, list=FALSE)
train <- integrated[train_idx,]
temp <- integrated[-train_idx,]

val_idx <- createDataPartition(temp$progression_group, p=0.50, list=FALSE)
val <- temp[val_idx,]
test <- temp[-val_idx,]

cat(sprintf("  Train: %d (Early: %d, Late: %d)\n", 
            nrow(train),
            sum(train$progression_group=="early"),
            sum(train$progression_group=="late")))
cat(sprintf("  Val:   %d (Early: %d, Late: %d)\n", 
            nrow(val),
            sum(val$progression_group=="early"),
            sum(val$progression_group=="late")))
cat(sprintf("  Test:  %d (Early: %d, Late: %d)\n\n", 
            nrow(test),
            sum(test$progression_group=="early"),
            sum(test$progression_group=="late")))

# Save
dir.create("data/integrated", showWarnings=FALSE, recursive=TRUE)
dir.create("data/splits", showWarnings=FALSE, recursive=TRUE)

saveRDS(integrated, "data/integrated/multi_modal_integrated.rds")
saveRDS(train, "data/splits/train_data.rds")
saveRDS(val, "data/splits/val_data.rds")
saveRDS(test, "data/splits/test_data.rds")

cat("=== MULTI-MODAL INTEGRATION COMPLETE ===\n\n")

# Feature summary
all_features <- colnames(integrated)[!colnames(integrated) %in% c("patient_barcode", "progression_group")]

cat("Feature counts by type:\n")
cat(sprintf("  Clinical: %d\n", sum(!str_detect(all_features, "PATHWAY_|miRNA_|INTERACT_"))))
cat(sprintf("  Pathway: %d\n", sum(str_detect(all_features, "^PATHWAY_"))))
cat(sprintf("  miRNA: %d\n", sum(str_detect(all_features, "^miRNA_"))))
cat(sprintf("  Interaction: %d\n", sum(str_detect(all_features, "^INTERACT_"))))
cat(sprintf("  Total: %d features\n\n", length(all_features)))

cat("Files saved:\n")
cat("  - data/integrated/multi_modal_integrated.rds\n")
cat("  - data/splits/train_data.rds\n")
cat("  - data/splits/val_data.rds\n")
cat("  - data/splits/test_data.rds\n\n")

cat("Next: Rscript scripts/07_eda_and_de_analysis.R\n\n")