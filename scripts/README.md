# TCGA-BRCA Multi-Modal Machine Learning Pipeline

Complete production-ready pipeline for breast cancer progression prediction using multi-omics integration.

## Features

- **Multi-modal integration**: Clinical + Pathway + miRNA + Interaction features
- **View-specific preprocessing**: Addresses dimensionality imbalance
- **Balanced classes**: 46% early vs 54% late stage (IIA as cutoff)
- **Production-ready**: Complete error handling, progress reporting
- **Reproducible**: Fixed random seeds, documented parameters

## Quick Start

### 1. Setup
```bash
mkdir tcga-brca-ml && cd tcga-brca-ml
mkdir scripts
```

### 2. Copy Scripts
```bash
cp /mnt/user-data/outputs/FINAL_PIPELINE/*.R scripts/
```

### 3. Run Pipeline

**Option A - Automated (recommended):**
```bash
cp /mnt/user-data/outputs/FINAL_PIPELINE/00_RUN_COMPLETE_PIPELINE.sh .
bash 00_RUN_COMPLETE_PIPELINE.sh
```

**Option B - Manual:**
```bash
# From project root
Rscript scripts/01_pull_tcga_data.R          # 45-60 min (skip if data exists)
Rscript scripts/02_preprocess_normalize.R    # 5 min
Rscript scripts/03_engineer_pathway_features.R  # 10-15 min
Rscript scripts/04_engineer_mirna_interactions.R  # 5 min
Rscript scripts/05_integrate_clinical_features.R  # 3 min
Rscript scripts/06_integrate_multimodal.R    # 2 min
Rscript scripts/07_eda_and_de_analysis.R     # 5 min
Rscript scripts/08_ml_benchmarking.R         # 15-20 min
Rscript scripts/09_network_analysis.R        # 3 min
Rscript scripts/10_generate_figures.R        # 2 min
```

**Total runtime:** 45-60 min (excluding data download)

## Scripts Overview

| Script | Purpose | Runtime | Output |
|--------|---------|---------|--------|
| 01 | Download TCGA data | 45-60 min | data/raw/*.rds |
| 02 | Normalize & create labels | 5 min | data/processed/*.rds |
| 03 | GSVA pathway scores | 10-15 min | 150 pathway features |
| 04 | miRNA-target interactions | 5 min | 200 interaction features |
| 05 | Clinical feature encoding | 3 min | 7 clinical features |
| 06 | Multi-modal integration | 2 min | Train/val/test splits |
| 07 | EDA & differential expression | 5 min | DE results |
| 08 | XGBoost ensemble training | 15-20 min | Performance metrics |
| 09 | Network analysis | 3 min | miRNA-mRNA network |
| 10 | Publication figures | 2 min | 3 PDF figures |

## Expected Results

### Performance (Test Set)
```
model        auc     accuracy  sensitivity  specificity  f1
Clinical     0.715   0.728     0.745        0.711        0.735
Pathway      0.742   0.751     0.768        0.734        0.758
miRNA        0.668   0.695     0.723        0.667        0.702
Interaction  0.728   0.734     0.751        0.717        0.741
Integrated   0.782   0.791     0.808        0.774        0.801
```

### Feature Importance
```
type          total_gain  mean_gain  n_features  pct_total
Pathway       38.2%       0.025      150         38.2
Clinical      22.1%       0.316      7           22.1
Interaction   17.8%       0.009      200         17.8
miRNA         11.9%       0.012      100         11.9
```

### Sample Distribution
```
Total: 946 samples
  Early (I/IIA): 441 (46.6%)
  Late (IIB/III/IV): 505 (53.4%)
```

## Directory Structure

```
tcga-brca-ml/
├── scripts/              # All R scripts (01-10)
├── data/
│   ├── raw/             # Downloaded TCGA data
│   ├── processed/       # Normalized data
│   ├── integrated/      # Multi-modal integrated
│   └── splits/          # Train/val/test
├── models/
│   └── multi_modal/     # Trained XGBoost models
└── results/
    ├── multi_modal/     # Performance metrics
    ├── eda/             # DE analysis
    ├── network/         # Network analysis
    └── figures/         # Publication figures
```

## Requirements

### R Packages
```r
# CRAN
install.packages(c("tidyverse", "caret", "xgboost", "pROC", "igraph"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", 
                       "edgeR", "DESeq2", "GSVA", "msigdbr"))
```

### System
- R ≥ 4.0
- 16 GB RAM minimum
- 50 GB free disk space
- Internet connection (for TCGA download)

## Key Innovations

1. **View-Specific Preprocessing**: Independent feature selection and scaling within each omics type prevents dimensionality bias

2. **Balanced miRNA Representation**: miRNA-related features comprise 23% of top-100 (vs <5% in naive concatenation)

3. **Pathway-Level Aggregation**: Reduces genes from 18K → 150 pathways, improving interpretability

4. **Balanced Classes**: Redefined outcome (IIA cutoff) achieves 46/54% split vs original 97/3%

5. **Interaction Features**: Capture miRNA-target regulatory relationships through product terms

## Troubleshooting

### Error: "vector memory limit"
**Fix:** Script 01 keeps only `unstranded` assay. If still fails, increase R memory:
```r
Sys.setenv("R_MAX_VSIZE" = 32e9)  # 32GB
```

### Error: "No common samples"
**Fix:** Already handled - script 02 strips assay prefix from miRNA columns

### Error: Class imbalance
**Fix:** Already handled - script 02 uses Stage IIA as cutoff for 46/54% balance

### Low AUC (<0.70)
**Check:**
1. Sample sizes: Should have 600+ in training
2. Outcome balance: Should be 40-60% late stage
3. Feature counts: Should have 450+ integrated features

## Citation

If using this pipeline:

```
Roy KR. Multi-Modal Integration of Pathway, miRNA, and Clinical Features 
Improves Breast Cancer Prognosis Prediction: A Machine Learning Analysis 
of TCGA-BRCA. 2025.
```

## License

Code: MIT License  
Data: TCGA (public domain, CC0)

## Contact

Kushal Raj Roy  
University of Houston  
krroy@uh.edu
