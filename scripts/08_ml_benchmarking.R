#!/usr/bin/env Rscript
################################################################################
# STEP 08: ML BENCHMARKING
################################################################################

pkgs <- c("xgboost", "caret", "pROC")
for(p in pkgs) {
  if(!require(p, quietly=TRUE, character.only=TRUE)) {
    install.packages(p, repos="https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(xgboost)
  library(caret)
  library(pROC)
})

cat("\n=== ML BENCHMARKING ===\n\n")

# Load splits
train <- readRDS("data/splits/train_data.rds")
val <- readRDS("data/splits/val_data.rds")
test <- readRDS("data/splits/test_data.rds")

# Separate features by type
all_feat <- colnames(train)[!colnames(train) %in% c("patient_barcode", "progression_group")]
clin_feat <- all_feat[!str_detect(all_feat, "PATHWAY_|miRNA_|INTERACT_")]
path_feat <- all_feat[str_detect(all_feat, "^PATHWAY_")]
mirna_feat <- all_feat[str_detect(all_feat, "^miRNA_")]
inter_feat <- all_feat[str_detect(all_feat, "^INTERACT_")]

cat(sprintf("Features:\n"))
cat(sprintf("  Clinical: %d\n", length(clin_feat)))
cat(sprintf("  Pathway: %d\n", length(path_feat)))
cat(sprintf("  miRNA: %d\n", length(mirna_feat)))
cat(sprintf("  Interaction: %d\n\n", length(inter_feat)))

# Outcomes
y_train <- as.numeric(train$progression_group == "late")
y_val <- as.numeric(val$progression_group == "late")
y_test <- as.numeric(test$progression_group == "late")

# Training function
train_xgb <- function(X_tr, y_tr, X_val, y_val, name) {
  cat(sprintf("Training %s model...\n", name))
  
  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    eta = 0.05,
    max_depth = 6,
    subsample = 0.8,
    colsample_bytree = 0.8
  )
  
  model <- xgb.train(
    params = params,
    data = xgb.DMatrix(as.matrix(X_tr), label=y_tr),
    nrounds = 500,
    watchlist = list(val = xgb.DMatrix(as.matrix(X_val), label=y_val)),
    early_stopping_rounds = 20,
    verbose = 0
  )
  
  return(model)
}

# Train models
models <- list(
  clinical = train_xgb(train[,clin_feat], y_train, val[,clin_feat], y_val, "Clinical"),
  pathway = train_xgb(train[,path_feat], y_train, val[,path_feat], y_val, "Pathway"),
  mirna = train_xgb(train[,mirna_feat], y_train, val[,mirna_feat], y_val, "miRNA"),
  interact = train_xgb(train[,inter_feat], y_train, val[,inter_feat], y_val, "Interaction"),
  integrated = train_xgb(train[,all_feat], y_train, val[,all_feat], y_val, "Integrated")
)

cat("\n")

# Evaluation function
evaluate <- function(model, X, y, name) {
  pred <- predict(model, xgb.DMatrix(as.matrix(X)))
  auc_val <- as.numeric(auc(roc(y, pred, quiet=TRUE)))
  pred_class <- ifelse(pred > 0.5, 1, 0)
  
  cm <- confusionMatrix(
    factor(pred_class, levels=c(0,1)),
    factor(y, levels=c(0,1))
  )
  
  data.frame(
    model = name,
    auc = auc_val,
    accuracy = mean(pred_class == y),
    sensitivity = as.numeric(cm$byClass["Sensitivity"]),
    specificity = as.numeric(cm$byClass["Specificity"]),
    f1 = as.numeric(cm$byClass["F1"])
  )
}

# Evaluate on test set
results <- bind_rows(
  evaluate(models$clinical, test[,clin_feat], y_test, "Clinical"),
  evaluate(models$pathway, test[,path_feat], y_test, "Pathway"),
  evaluate(models$mirna, test[,mirna_feat], y_test, "miRNA"),
  evaluate(models$interact, test[,inter_feat], y_test, "Interaction"),
  evaluate(models$integrated, test[,all_feat], y_test, "Integrated")
)

# Save results
dir.create("results/multi_modal", showWarnings=FALSE, recursive=TRUE)
dir.create("models/multi_modal", showWarnings=FALSE, recursive=TRUE)

write_csv(results, "results/multi_modal/performance_summary.csv")
xgb.save(models$integrated, "models/multi_modal/model_integrated.xgb")

cat("TEST SET PERFORMANCE:\n")
print(results, digits=3)
cat("\n")

# Feature importance
importance <- xgb.importance(model=models$integrated)
importance$type <- case_when(
  str_detect(importance$Feature, "^PATHWAY_") ~ "Pathway",
  str_detect(importance$Feature, "^miRNA_") ~ "miRNA",
  str_detect(importance$Feature, "^INTERACT_") ~ "Interaction",
  TRUE ~ "Clinical"
)

write_csv(importance, "results/multi_modal/feature_importance.csv")

importance_by_type <- importance %>%
  group_by(type) %>%
  summarise(
    total_gain = sum(Gain),
    mean_gain = mean(Gain),
    n_features = n(),
    pct_total = 100 * sum(Gain) / sum(importance$Gain)
  ) %>%
  arrange(desc(total_gain))

write_csv(importance_by_type, "results/multi_modal/importance_by_type.csv")

cat("Importance by feature type:\n")
print(importance_by_type, digits=3)
cat("\n")

cat("✓ Complete\n\n")
