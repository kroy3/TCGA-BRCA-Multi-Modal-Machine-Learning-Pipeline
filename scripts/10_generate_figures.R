#!/usr/bin/env Rscript
################################################################################
# STEP 10: GENERATE PUBLICATION FIGURES
################################################################################

if(!require(xgboost, quietly=TRUE)) {
  install.packages("xgboost", repos="https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(xgboost)
})

cat("\n=== GENERATING FIGURES ===\n\n")

# Load data
clinical <- readRDS("data/processed/clinical_integrated/clinical_ml_features.rds")
perf <- read_csv("results/multi_modal/performance_summary.csv", show_col_types=FALSE)

dir.create("results/figures", showWarnings=FALSE, recursive=TRUE)

## Figure 1: Cohort Distribution
cat("Creating Figure 1: Cohort distribution...\n")

cohort <- data.frame(
  stage = c("Early\n(I/IIA)", "Late\n(IIB/III/IV)"),
  n = c(
    sum(clinical$progression_group=="early", na.rm=TRUE),
    sum(clinical$progression_group=="late", na.rm=TRUE)
  )
)

pdf("results/figures/figure1_cohort.pdf", width=7, height=5)
par(mar=c(5,5,3,2))
bp <- barplot(
  cohort$n,
  names.arg = cohort$stage,
  main = "TCGA-BRCA Cohort Distribution",
  ylab = "Number of Samples",
  col = c("steelblue3", "coral3"),
  ylim = c(0, max(cohort$n)*1.2),
  cex.names = 1.2,
  cex.axis = 1.1,
  cex.lab = 1.2
)
text(bp, cohort$n, labels=cohort$n, pos=3, cex=1.3, font=2)
dev.off()

## Figure 2: Model Performance
cat("Creating Figure 2: Model performance...\n")

pdf("results/figures/figure2_performance.pdf", width=9, height=6)
par(mar=c(9,5,3,2))
bp <- barplot(
  perf$auc,
  names.arg = perf$model,
  main = "Model Performance Comparison (Test Set)",
  ylab = "AUC",
  ylim = c(0, 1),
  col = "steelblue3",
  las = 2,
  cex.names = 1.1,
  cex.axis = 1.1,
  cex.lab = 1.2
)
abline(h=0.75, lty=2, col="red3", lwd=2)
text(bp, perf$auc, labels=sprintf("%.3f", perf$auc), pos=3, cex=1.0)
legend("bottomright", legend="AUC = 0.75", lty=2, col="red3", lwd=2, cex=1.1)
dev.off()

## Figure 3: Feature Importance
cat("Creating Figure 3: Feature importance...\n")

model <- xgb.load("models/multi_modal/model_integrated.xgb")
importance <- xgb.importance(model=model)

top30 <- head(importance, 30)
top30$type <- case_when(
  str_detect(top30$Feature, "^PATHWAY_") ~ "Pathway",
  str_detect(top30$Feature, "^miRNA_") ~ "miRNA",
  str_detect(top30$Feature, "^INTERACT_") ~ "Interaction",
  TRUE ~ "Clinical"
)

# Clean feature names
top30$Feature_clean <- str_replace(top30$Feature, "^(PATHWAY_|miRNA_|INTERACT_)", "")
top30$Feature_clean <- str_trunc(top30$Feature_clean, 40)

pdf("results/figures/figure3_importance.pdf", width=11, height=10)
par(mar=c(5,15,3,2))

colors <- c(
  "Clinical" = "darkgreen",
  "Pathway" = "steelblue3",
  "miRNA" = "coral3",
  "Interaction" = "purple3"
)

bp <- barplot(
  rev(top30$Gain),
  names.arg = rev(top30$Feature_clean),
  horiz = TRUE,
  las = 1,
  col = colors[rev(top30$type)],
  main = "Top 30 Features by Importance",
  xlab = "Gain",
  cex.names = 0.7,
  cex.axis = 1.0,
  cex.lab = 1.1
)

legend("bottomright", 
       legend = names(colors),
       fill = colors,
       cex = 1.0,
       title = "Feature Type")
dev.off()

cat("\n✓ Generated 3 figures\n\n")

cat("=== PIPELINE COMPLETE ===\n\n")
cat("Results Summary:\n")
cat("  Performance: results/multi_modal/performance_summary.csv\n")
cat("  Importance: results/multi_modal/importance_by_type.csv\n")
cat("  Figures: results/figures/\n")
cat("  Models: models/multi_modal/\n\n")

cat("View results:\n")
cat("  cat results/multi_modal/performance_summary.csv\n")
cat("  cat results/multi_modal/importance_by_type.csv\n\n")
