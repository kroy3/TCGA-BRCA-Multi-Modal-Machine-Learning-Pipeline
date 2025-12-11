#!/usr/bin/env Rscript
################################################################################
# STEP 04: miRNA-TARGET INTERACTION FEATURES
# Fixed: Recover rownames and convert character to numeric
################################################################################

suppressPackageStartupMessages(library(tidyverse))

cat("\n=== miRNA-TARGET INTERACTIONS ===\n\n")

mirna_raw <- readRDS("data/processed/mirna_log2cpm.rds")
mrna_raw <- readRDS("data/processed/mrna_log2cpm.rds")

cat(sprintf("Loaded: %d miRNAs x %d samples, %d mRNAs x %d samples\n", 
            nrow(mirna_raw), ncol(mirna_raw), nrow(mrna_raw), ncol(mrna_raw)))

# FIX 1: Convert miRNA from character to numeric
cat("\nFixing miRNA data...\n")
cat(sprintf("  miRNA class: %s\n", class(mirna_raw)[1]))

if(is.character(mirna_raw[1,1])) {
  cat("  Converting from character to numeric\n")
  mirna <- matrix(as.numeric(mirna_raw), nrow=nrow(mirna_raw), ncol=ncol(mirna_raw))
  colnames(mirna) <- colnames(mirna_raw)
} else {
  mirna <- mirna_raw
}

# FIX 2: Recover or create miRNA rownames
if(is.null(rownames(mirna_raw))) {
  cat("  No rownames found - attempting recovery\n")
  
  # Try loading original file to get feature names
  if(file.exists("data/raw/mirna_raw.rds")) {
    cat("  Loading original miRNA data for feature names\n")
    mirna_original <- readRDS("data/raw/mirna_raw.rds")
    
    if(is(mirna_original, "SummarizedExperiment")) {
      library(SummarizedExperiment)
      mirna_features <- rownames(mirna_original)
    } else {
      mirna_features <- rownames(mirna_original)
    }
    
    if(!is.null(mirna_features) && length(mirna_features) == nrow(mirna)) {
      rownames(mirna) <- mirna_features
      cat(sprintf("  Recovered %d feature names\n", length(mirna_features)))
    } else {
      # Fallback: create generic names
      rownames(mirna) <- paste0("miRNA_", 1:nrow(mirna))
      cat(sprintf("  Created generic names: miRNA_1 to miRNA_%d\n", nrow(mirna)))
    }
  } else {
    # Create generic names
    rownames(mirna) <- paste0("miRNA_", 1:nrow(mirna))
    cat(sprintf("  Created generic names: miRNA_1 to miRNA_%d\n", nrow(mirna)))
  }
} else {
  rownames(mirna) <- rownames(mirna_raw)
}

# Ensure mRNA is numeric
mrna <- apply(mrna_raw, 2, as.numeric)
rownames(mrna) <- rownames(mrna_raw)
colnames(mrna) <- colnames(mrna_raw)

cat(sprintf("\nAfter fixing:\n"))
cat(sprintf("  miRNA: %d features with rownames\n", sum(!is.na(rownames(mirna)))))
cat(sprintf("  mRNA: %d features with rownames\n", sum(!is.na(rownames(mrna)))))
cat(sprintf("  First miRNA: %s\n", rownames(mirna)[1]))
cat(sprintf("  First mRNA: %s\n\n", rownames(mrna)[1]))

# Extract barcodes and match samples
extract_bc <- function(x) str_sub(as.character(x), 1, 16)

mirna_bc <- extract_bc(colnames(mirna))
mrna_bc <- extract_bc(colnames(mrna))
common_bc <- intersect(mirna_bc, mrna_bc)

cat(sprintf("Common samples: %d\n\n", length(common_bc)))

# Match
mirna_idx <- match(common_bc, mirna_bc)
mrna_idx <- match(common_bc, mrna_bc)

mirna_matched <- mirna[, mirna_idx, drop=FALSE]
mrna_matched <- mrna[, mrna_idx, drop=FALSE]

colnames(mirna_matched) <- common_bc
colnames(mrna_matched) <- common_bc

cat("Normalizing miRNA...\n")

# Log transform if needed
mirna_max <- max(mirna_matched, na.rm=TRUE)
cat(sprintf("  Max value: %.2f\n", mirna_max))

if(mirna_max > 100) {
  cat("  Applying log2(x + 1)\n")
  mirna_matched <- log2(mirna_matched + 1)
}

# Standardize
mirna_scaled <- t(scale(t(mirna_matched)))
rownames(mirna_scaled) <- rownames(mirna_matched)
colnames(mirna_scaled) <- common_bc

cat(sprintf("  Final range: %.2f to %.2f\n\n", 
            min(mirna_scaled, na.rm=TRUE), 
            max(mirna_scaled, na.rm=TRUE)))

cat("Selecting variable features...\n")

# Compute SD
mirna_sd <- apply(mirna_scaled, 1, sd, na.rm=TRUE)
mrna_sd <- apply(mrna_matched, 1, sd, na.rm=TRUE)

mirna_sd[is.na(mirna_sd)] <- 0
mrna_sd[is.na(mrna_sd)] <- 0

# Select top
top_mirna <- names(sort(mirna_sd, decreasing=TRUE)[1:min(50, sum(mirna_sd > 0))])
top_mrna <- names(sort(mrna_sd, decreasing=TRUE)[1:min(500, sum(mrna_sd > 0))])

cat(sprintf("  Selected: %d miRNAs x %d mRNAs\n\n", length(top_mirna), length(top_mrna)))

if(length(top_mirna) == 0) {
  cat("ERROR: No variable miRNAs\n")
  quit(status=1)
}

cat("Computing interactions...\n")

# Subsets
mirna_subset <- mirna_scaled[top_mirna, , drop=FALSE]
mrna_subset <- mrna_matched[top_mrna, , drop=FALSE]

# Compute
interactions <- matrix(0, nrow=200, ncol=length(common_bc))
colnames(interactions) <- common_bc
idx <- 1

for(i in 1:nrow(mirna_subset)) {
  if(idx > 200) break
  
  mirna_expr <- as.numeric(mirna_subset[i, ])
  
  cors <- apply(mrna_subset, 1, function(x) {
    cor(mirna_expr, as.numeric(x), method="spearman", use="complete.obs")
  })
  
  strong_neg <- which(cors < -0.3 & !is.na(cors))
  
  if(length(strong_neg) > 0) {
    for(j in strong_neg[1:min(4, length(strong_neg))]) {
      if(idx > 200) break
      mrna_expr <- as.numeric(mrna_subset[j, ])
      interactions[idx,] <- scale(mirna_expr)[,1] * scale(mrna_expr)[,1]
      idx <- idx + 1
    }
  }
  
  if(i %% 10 == 0) cat(sprintf("  Progress: %d/%d\n", i, nrow(mirna_subset)))
}

cat(sprintf("\nCreated %d interactions\n\n", idx - 1))

# Finalize
if(idx > 1) {
  interactions <- interactions[1:(idx-1), , drop=FALSE]
} else {
  n <- min(50, nrow(mirna_subset))
  interactions <- matrix(0, nrow=n, ncol=length(common_bc))
  colnames(interactions) <- common_bc
  for(i in 1:n) {
    interactions[i,] <- scale(as.numeric(mirna_subset[i,]))[,1] * 
                        scale(as.numeric(mrna_subset[i,]))[,1]
  }
}

rownames(interactions) <- paste0("INTERACT_", 1:nrow(interactions))

# Save
dir.create("data/processed/mirna_target_features", showWarnings=FALSE, recursive=TRUE)
saveRDS(interactions, "data/processed/mirna_target_features/mirna_target_interactions.rds")

cat("=== COMPLETE ===\n\n")
cat(sprintf("Final: %d interactions x %d samples\n\n", 
            nrow(interactions), ncol(interactions)))

cat("Next: Rscript scripts/05_integrate_clinical_features.R\n\n")