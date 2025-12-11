#!/bin/bash
################################################################################
# MASTER PIPELINE RUNNER
# Executes all 10 scripts in sequence
################################################################################

echo "=========================================="
echo "  TCGA-BRCA ML PIPELINE"
echo "=========================================="
echo ""

# Check we're in project root
if [ ! -d "scripts" ]; then
    echo "ERROR: Run from project root (not scripts/)"
    exit 1
fi

START_TIME=$(date +%s)

# Function to run script
run_script() {
    local script=$1
    local name=$2
    
    echo "--------------------------------------"
    echo "Running: $name"
    echo "--------------------------------------"
    
    if Rscript scripts/$script; then
        echo "✓ $name complete"
        echo ""
    else
        echo "✗ $name FAILED"
        exit 1
    fi
}

# Run pipeline (skip 01 if data exists)
if [ -f "data/raw/mirna_raw.rds" ]; then
    echo "⊘ Skipping data download (files exist)"
    echo ""
else
    run_script "01_pull_tcga_data.R" "Data Download"
fi

run_script "02_preprocess_normalize.R" "Preprocessing"
run_script "03_engineer_pathway_features.R" "Pathway Engineering"
run_script "04_engineer_mirna_interactions.R" "miRNA Interactions"
run_script "05_integrate_clinical_features.R" "Clinical Integration"
run_script "06_integrate_multimodal.R" "Multi-Modal Integration"
run_script "07_eda_and_de_analysis.R" "EDA & DE Analysis"
run_script "08_ml_benchmarking.R" "ML Benchmarking"
run_script "09_network_analysis.R" "Network Analysis"
run_script "10_generate_figures.R" "Figure Generation"

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo "=========================================="
echo "  PIPELINE COMPLETE"
echo "=========================================="
echo ""
echo "Total runtime: $((RUNTIME / 60)) minutes"
echo ""
echo "Results:"
echo "  cat results/multi_modal/performance_summary.csv"
echo "  cat results/multi_modal/importance_by_type.csv"
echo ""
echo "Figures:"
echo "  open results/figures/*.pdf"
echo ""
