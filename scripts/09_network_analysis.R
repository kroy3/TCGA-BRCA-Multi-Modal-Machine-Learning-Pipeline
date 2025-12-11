#!/usr/bin/env Rscript
################################################################################
# STEP 09: NETWORK ANALYSIS (final, robust with diagnostics)
################################################################################

if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
})

cat("\n=== NETWORK ANALYSIS ===\n\n")

# -------------------------------------------------------------------------
# 1. Load and coerce to numeric
# -------------------------------------------------------------------------
mirna_raw <- readRDS("data/processed/mirna_log2cpm.rds")
mrna_raw  <- readRDS("data/processed/mrna_log2cpm.rds")

mirna <- as.matrix(mirna_raw)
mrna  <- as.matrix(mrna_raw)

mode(mirna) <- "numeric"
mode(mrna)  <- "numeric"

# Give miRNA rows explicit names if they don't have any
if (is.null(rownames(mirna))) {
  rownames(mirna) <- paste0("miRNA_", seq_len(nrow(mirna)))
}

if (any(!is.finite(mirna))) {
  warning("Non-finite values in miRNA matrix; setting to NA")
  mirna[!is.finite(mirna)] <- NA
}
if (any(!is.finite(mrna))) {
  warning("Non-finite values in mRNA matrix; setting to NA")
  mrna[!is.finite(mrna)] <- NA
}

# -------------------------------------------------------------------------
# 2. Align samples between miRNA and mRNA
# -------------------------------------------------------------------------
common <- intersect(colnames(mirna), colnames(mrna))

if (length(common) == 0) {
  cat("No exact matching column names; trying TCGA-style 12-char patient barcodes...\n")
  extract_bc <- function(x) stringr::str_sub(as.character(x), 1, 12)

  mirna_bc <- extract_bc(colnames(mirna))
  mrna_bc  <- extract_bc(colnames(mrna))

  common_bc <- intersect(mirna_bc, mrna_bc)

  if (length(common_bc) == 0) {
    stop("No overlapping samples or patient barcodes between miRNA and mRNA matrices. Check column names.")
  }

  common_bc <- sort(common_bc)

  mirna_idx <- match(common_bc, mirna_bc)
  mrna_idx  <- match(common_bc, mrna_bc)

  mirna <- mirna[, mirna_idx, drop = FALSE]
  mrna  <- mrna[, mrna_idx,  drop = FALSE]

  colnames(mirna) <- common_bc
  colnames(mrna)  <- common_bc
  common <- colnames(mirna)
} else {
  mirna <- mirna[, common, drop = FALSE]
  mrna  <- mrna[, common, drop = FALSE]
}

n_samples <- length(common)
cat(sprintf("Computing correlations on %d samples\n", n_samples))
if (n_samples == 0) {
  stop("After alignment, there are still 0 common samples between miRNA and mRNA.")
}

# -------------------------------------------------------------------------
# 3. Remove non-variable features
# -------------------------------------------------------------------------
var_mirna <- apply(mirna, 1, var, na.rm = TRUE)
var_mrna  <- apply(mrna, 1, var, na.rm = TRUE)

mirna <- mirna[is.finite(var_mirna) & var_mirna > 0, , drop = FALSE]
mrna  <- mrna[is.finite(var_mrna)  & var_mrna  > 0, , drop = FALSE]

cat(sprintf("  miRNA features after variance filter: %d\n", nrow(mirna)))
cat(sprintf("  mRNA  features after variance filter: %d\n", nrow(mrna)))

if (nrow(mirna) == 0 || nrow(mrna) == 0) {
  stop("No variable miRNA or mRNA features remain after variance filtering.")
}

# Re-slice variance vectors to current rows and keep names
var_mirna <- var_mirna[rownames(mirna)]
var_mrna  <- var_mrna[rownames(mrna)]

# -------------------------------------------------------------------------
# 4. Choose top variable miRNAs and mRNAs
# -------------------------------------------------------------------------
n_top_mirna <- min(20, nrow(mirna))
n_top_mrna  <- min(500, nrow(mrna))

if (n_top_mirna == 0) stop("No miRNA features available for network analysis after filtering.")
if (n_top_mrna  == 0) stop("No mRNA features available for network analysis after filtering.")

top_mirna <- names(sort(var_mirna, decreasing = TRUE))[seq_len(n_top_mirna)]
top_mrna  <- names(sort(var_mrna,  decreasing = TRUE))[seq_len(n_top_mrna)]

cat(sprintf("  Using %d miRNAs and %d mRNAs for network\n",
            length(top_mirna), length(top_mrna)))

# -------------------------------------------------------------------------
# 5. Correlation matrix
# -------------------------------------------------------------------------
cat("Computing miRNA–mRNA correlation matrix...\n")

cor_mat <- cor(
  t(mirna[top_mirna, , drop = FALSE]),
  t(mrna[top_mrna,  , drop = FALSE]),
  use = "pairwise.complete.obs"
)

# Diagnostics on correlation distribution
cor_vals   <- as.vector(cor_mat)
cor_range  <- range(cor_vals, na.rm = TRUE)
max_abscor <- max(abs(cor_vals), na.rm = TRUE)

cat(sprintf("Correlation range: %.3f to %.3f; max |cor| = %.3f\n",
            cor_range[1], cor_range[2], max_abscor))

# -------------------------------------------------------------------------
# 6. Extract strong edges (with fallback)
# -------------------------------------------------------------------------
cor_cutoff <- 0.30  # more lenient than 0.5

strong_idx <- which(abs(cor_mat) > cor_cutoff & !is.na(cor_mat), arr.ind = TRUE)

if (nrow(strong_idx) == 0) {
  cat(sprintf("No pairs with |cor| > %.2f. Taking top 200 pairs by |cor| instead.\n",
              cor_cutoff))

  # Order all pairs by absolute correlation
  ord   <- order(abs(cor_vals), decreasing = TRUE)
  top_k <- min(200, length(ord))

  # Convert vector indices back to matrix row/col indices
  strong_idx <- arrayInd(ord[seq_len(top_k)], dim(cor_mat))
}

edges <- tibble(
  mirna       = rownames(cor_mat)[strong_idx[, 1]],
  mrna        = colnames(cor_mat)[strong_idx[, 2]],
  correlation = cor_mat[strong_idx]
)

# -------------------------------------------------------------------------
# 7. Build graph and compute metrics
# -------------------------------------------------------------------------
cat(sprintf("Found %d selected edges. Building graph...\n", nrow(edges)))

g <- graph_from_data_frame(edges, directed = FALSE)

metrics <- tibble(
  feature     = V(g)$name,
  type        = case_when(
    V(g)$name %in% top_mirna ~ "miRNA",
    V(g)$name %in% top_mrna  ~ "mRNA",
    TRUE                     ~ "unknown"
  ),
  degree      = degree(g),
  betweenness = betweenness(g),
  closeness   = closeness(g)
) %>%
  arrange(desc(degree))

# -------------------------------------------------------------------------
# 8. Save results
# -------------------------------------------------------------------------
dir.create("results/network", showWarnings = FALSE, recursive = TRUE)

readr::write_csv(edges,   "results/network/network_edges.csv")
readr::write_csv(metrics, "results/network/network_metrics.csv")

saveRDS(g, "results/network/network_graph.rds")

cat(sprintf("\nNetwork: %d edges, %d nodes\n\n", nrow(edges), vcount(g)))

cat("Top 10 hub features:\n")
print(head(metrics, 10))
cat("\n")

cat("✓ Complete\n\n")
