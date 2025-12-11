#!/usr/bin/env Rscript
################################################################################
# STEP 03: PATHWAY FEATURE ENGINEERING
# Fixed: Use correct msigdbr 10.0+ column names
################################################################################

if(!require(GSVA, quietly=TRUE)) {
  if(!require(BiocManager, quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install("GSVA")
}

if(!require(BiocParallel, quietly=TRUE)) {
  if(!require(BiocManager, quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install("BiocParallel")
}

if(!require(org.Hs.eg.db, quietly=TRUE)) {
  if(!require(BiocManager, quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
}

if(!require(msigdbr, quietly=TRUE)) install.packages("msigdbr")

suppressPackageStartupMessages({
  library(tidyverse)
  library(GSVA)
  library(BiocParallel)
  library(msigdbr)
  library(org.Hs.eg.db)
})

cat("\n=== PATHWAY ENGINEERING ===\n\n")

mrna <- readRDS("data/processed/mrna_log2cpm.rds")
cat(sprintf("mRNA: %d genes x %d samples\n", nrow(mrna), ncol(mrna)))

# Convert Ensembl IDs to gene symbols
cat("Converting gene IDs...\n")
gene_ids <- rownames(mrna)

if(any(grepl("^ENSG", gene_ids))) {
  cat("  Detected Ensembl IDs, converting to symbols...\n")
  
  gene_ids_clean <- sub("\\..*", "", gene_ids)
  
  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_ids_clean,
    columns = c("SYMBOL"),
    keytype = "ENSEMBL"
  )
  
  gene_map_df <- as.data.frame(gene_map)
  gene_map_df <- gene_map_df[!is.na(gene_map_df$SYMBOL), ]
  gene_map_df <- gene_map_df[!duplicated(gene_map_df$ENSEMBL), ]
  
  idx <- match(gene_ids_clean, gene_map_df$ENSEMBL)
  new_names <- gene_map_df$SYMBOL[idx]
  
  keep <- !is.na(new_names)
  mrna_filtered <- mrna[keep, ]
  rownames(mrna_filtered) <- new_names[keep]
  mrna_filtered <- mrna_filtered[!duplicated(rownames(mrna_filtered)), ]
  
  cat(sprintf("  Converted: %d -> %d genes with symbols\n\n", nrow(mrna), nrow(mrna_filtered)))
  
} else {
  mrna_filtered <- mrna
}

# Load pathways
cat("Loading pathways...\n")

# 1. Hallmark (50 pathways)
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
cat(sprintf("  Hallmark: %d pathways\n", length(hallmark_sets)))

# 2. Reactome - check column names first
reactome_all <- msigdbr(species = "Homo sapiens", collection = "C2")

# Check available columns
cat("  Checking C2 columns: ", paste(names(reactome_all), collapse=", "), "\n")

# Filter for Reactome (use correct column name)
if("gs_subcollection" %in% names(reactome_all)) {
  reactome <- reactome_all[reactome_all$gs_subcollection == "CP:REACTOME", ]
} else if("gs_subcat" %in% names(reactome_all)) {
  reactome <- reactome_all[reactome_all$gs_subcat == "CP:REACTOME", ]
} else {
  # Fallback: filter by name prefix
  reactome <- reactome_all[grepl("^REACTOME_", reactome_all$gs_name), ]
}

reactome_sets <- split(reactome$gene_symbol, reactome$gs_name)
cat(sprintf("  Reactome: %d pathways\n", length(reactome_sets)))

# 3. GO Biological Process
gobp_all <- msigdbr(species = "Homo sapiens", collection = "C5")

if("gs_subcollection" %in% names(gobp_all)) {
  gobp <- gobp_all[gobp_all$gs_subcollection == "GO:BP", ]
} else if("gs_subcat" %in% names(gobp_all)) {
  gobp <- gobp_all[gobp_all$gs_subcat == "GO:BP", ]
} else {
  gobp <- gobp_all[grepl("^GOBP_", gobp_all$gs_name), ]
}

# Filter to reasonable sizes (20-500 genes)
gobp_sizes <- table(gobp$gs_name)
gobp_keep <- names(gobp_sizes)[gobp_sizes >= 20 & gobp_sizes <= 500]
gobp <- gobp[gobp$gs_name %in% gobp_keep, ]
gobp_sets <- split(gobp$gene_symbol, gobp$gs_name)
cat(sprintf("  GO:BP (filtered): %d pathways\n\n", length(gobp_sets)))

# Combine
all_gene_sets <- c(hallmark_sets, reactome_sets, gobp_sets)
cat(sprintf("Total gene sets: %d\n\n", length(all_gene_sets)))

# Compute GSVA
cat("Computing GSVA scores (15-25 min)...\n\n")

mrna_matrix <- as.matrix(mrna_filtered)

cat("  Processing all pathways...\n")
pathway_param <- gsvaParam(mrna_matrix, all_gene_sets, 
                          kcdf="Gaussian", minSize=15)
all_scores <- gsva(pathway_param, verbose=FALSE, BPPARAM=SerialParam())

cat(sprintf("\n  Computed: %d pathway scores\n", nrow(all_scores)))

# Select top 150 by variance
pathway_var <- apply(all_scores, 1, var)
top150 <- names(sort(pathway_var, decreasing=TRUE)[1:min(150, length(pathway_var))])
pathway_selected <- all_scores[top150,]

# Categorize
categorize_pathway <- function(name) {
  if(grepl("^HALLMARK", name)) return("Hallmark")
  if(grepl("^REACTOME", name)) return("Reactome")
  if(grepl("^GOBP", name)) return("GO:BP")
  return("Other")
}

# Save
dir.create("data/processed/pathway_features", showWarnings=FALSE, recursive=TRUE)
saveRDS(pathway_selected, "data/processed/pathway_features/gsva_pathway_scores.rds")

pathway_metadata <- data.frame(
  pathway = rownames(pathway_selected),
  variance = pathway_var[top150],
  mean_score = rowMeans(pathway_selected),
  source = sapply(rownames(pathway_selected), categorize_pathway)
)
pathway_metadata <- pathway_metadata[order(-pathway_metadata$variance), ]

write_csv(pathway_metadata, "data/processed/pathway_features/pathway_metadata.csv")

cat("\n=== PATHWAY ENGINEERING COMPLETE ===\n\n")
cat(sprintf("Selected: %d pathways x %d samples\n\n", 
            nrow(pathway_selected), ncol(pathway_selected)))

cat("Pathway sources:\n")
print(table(pathway_metadata$source))
cat("\n")

cat("Top 10 most variable pathways:\n")
top10 <- head(pathway_metadata, 10)
top10$pathway_short <- substr(top10$pathway, 1, 70)
print(top10[, c("pathway_short", "variance", "source")])
cat("\n")

# Check cancer-relevant
cancer_keywords <- c("CELL_CYCLE", "PROLIFERATION", "APOPTOSIS", "DNA_REPAIR", 
                     "PI3K", "MAPK", "P53", "MYC", "IMMUNE", "INFLAMMATION",
                     "E2F", "G2M")

cancer_related <- pathway_metadata[
  grepl(paste(cancer_keywords, collapse="|"), pathway_metadata$pathway), 
]

if(nrow(cancer_related) > 0) {
  cat(sprintf("Cancer-relevant pathways: %d\n", nrow(cancer_related)))
  cat("Examples:\n")
  print(head(cancer_related[, c("pathway", "source")], 10))
  cat("\n")
}

cat("Files saved:\n")
cat("  - data/processed/pathway_features/gsva_pathway_scores.rds\n")
cat("  - data/processed/pathway_features/pathway_metadata.csv\n\n")

cat("Next: Rscript scripts/04_engineer_mirna_interactions.R\n\n")