#=============================================================================
# 02_model_validation_cohort3.R
# Cohort 3 Disease Control Validation — HNSCC vs HC vs OPMDs
#=============================================================================
# Author:   Hexin Ma, Xinran Zhao
# Date:     2026-05
#
# Purpose:  Validates the 8-metabolite RF model on Cohort 3 (n=269):
#           (A) HNSCC vs HC — independent external validation
#           (B) HNSCC vs OPMDs — disease specificity
#           (C) OPMDs vs HC — selectivity check
#
# Input:    data/cohort1_8metabolites.csv
#           data/cohort3_8metabolites.csv
#           data/selected_features.csv
#
# Output:   output/tables/Table_cohort3_ROC.csv
#           output/tables/Table_cohort3_pred_scores.csv
#           output/figures/Fig_ROC_cohort3.pdf
#
# Usage:    Rscript R/02_model_validation_cohort3.R   (from project root)
# License:  GPL-3.0
#=============================================================================

source("R/00_setup.R")

library(randomForest)
library(pROC)
library(ROCR)

cat("=============================================================")
cat("02  Cohort 3 Disease Control Validation")
cat("=============================================================")

# ── 1. Load Data ─────────────────────────────────────────────────────────

cat("── Step 1: Loading data ──")

features <- read.csv("data/selected_features.csv")
metabolites <- features$Feature

cohort1 <- read.csv("data/cohort1_8metabolites.csv")
cohort3 <- read.csv("data/cohort3_8metabolites.csv")

cat(sprintf("  Cohort 1 (training): %d samples", nrow(cohort1)))
cat(sprintf("  Cohort 3 (disease control): %d samples", nrow(cohort3)))
cat("  Cohort 3 Group distribution:")
print(table(cohort3$Group))

# ── 2. Independent Z-Score Normalisation ──────────────────────────────────

cat("── Step 2: Independent Z-score normalisation ──")

# Cohort 1: Z-score using its own parameters
for (m in metabolites) {
  mu <- mean(cohort1[[m]], na.rm = TRUE)
  s  <- sd(cohort1[[m]], na.rm = TRUE)
  cohort1[[m]] <- (cohort1[[m]] - mu) / s
}

# Cohort 3: Z-score using its own parameters (all groups combined)
for (m in metabolites) {
  mu <- mean(cohort3[[m]], na.rm = TRUE)
  s  <- sd(cohort3[[m]], na.rm = TRUE)
  cohort3[[m]] <- (cohort3[[m]] - mu) / s
}

cat("  Cohort 1 post-Z mean (should be ~0):")
print(round(colMeans(cohort1[, metabolites]), 3))
cat("  Cohort 3 post-Z mean (should be ~0):")
print(round(colMeans(cohort3[, metabolites]), 3))

# ── 3. Train RF on Cohort 1 ──────────────────────────────────────────────

cat("── Step 3: Training RF on Cohort 1 (Z-scored) ──")

y_train <- factor(ifelse(cohort1$Group == "T", "HNSCC", "HC"),
                  levels = c("HC", "HNSCC"))

train_df <- cohort1[, metabolites, drop = FALSE]
train_df$Label <- y_train

set.seed(42)
rf_model <- randomForest(Label ~ ., data = train_df, ntree = 500, importance = TRUE)
cat("  RF trained (ntree=500)")

# ── 4. Predict All Cohort 3 Samples ──────────────────────────────────────

cat("── Step 4: Predicting Cohort 3 ──")

test_df <- cohort3[, metabolites, drop = FALSE]
pred_prob <- predict(rf_model, newdata = test_df, type = "prob")[, "HNSCC"]

cat(sprintf("  Predictions generated for %d samples", length(pred_prob)))

# ── 5. Analysis A: HNSCC vs HC ───────────────────────────────────────────

cat("── Step 5: HNSCC vs HC ──")

idx_A <- cohort3$Group %in% c("HNSCC", "HC")
y_A   <- ifelse(cohort3$Group[idx_A] == "HNSCC", 1, 0)
p_A   <- pred_prob[idx_A]

roc_A <- roc(y_A, p_A, levels = c(0, 1), direction = "<", quiet = TRUE)
ci_A  <- ci.auc(roc_A)
best_A <- coords(roc_A, x = "best", best.method = "youden",
                 ret = c("sensitivity", "specificity", "accuracy"))

cat(sprintf("  HNSCC vs HC (n=%d): AUC = %.4f (%.4f-%.4f)",
            sum(idx_A), ci_A[2], ci_A[1], ci_A[3]))
cat(sprintf("  Sensitivity = %.1f%%, Specificity = %.1f%%",
            best_A$sensitivity * 100, best_A$specificity * 100))

# ── 6. Analysis B: HNSCC vs OPMDs ────────────────────────────────────────

cat("── Step 6: HNSCC vs OPMDs ──")

idx_B <- cohort3$Group %in% c("HNSCC", "OPMDs")
y_B   <- ifelse(cohort3$Group[idx_B] == "HNSCC", 1, 0)
p_B   <- pred_prob[idx_B]

roc_B <- roc(y_B, p_B, levels = c(0, 1), direction = "<", quiet = TRUE)
ci_B  <- ci.auc(roc_B)

cat(sprintf("  HNSCC vs OPMDs (n=%d): AUC = %.4f (%.4f-%.4f)",
            sum(idx_B), ci_B[2], ci_B[1], ci_B[3]))

# ── 7. Analysis C: OPMDs vs HC ───────────────────────────────────────────

cat("── Step 7: OPMDs vs HC ──")

idx_C <- cohort3$Group %in% c("OPMDs", "HC")
y_C   <- ifelse(cohort3$Group[idx_C] == "OPMDs", 1, 0)
p_C   <- pred_prob[idx_C]

roc_C <- roc(y_C, p_C, levels = c(0, 1), direction = "<", quiet = TRUE)
ci_C  <- ci.auc(roc_C)

cat(sprintf("  OPMDs vs HC (n=%d): AUC = %.4f (%.4f-%.4f)",
            sum(idx_C), ci_C[2], ci_C[1], ci_C[3]))
cat("  (Low AUC expected — model is specific to HNSCC)")

# ── 8. Export Results ─────────────────────────────────────────────────────

cat("── Step 8: Exporting results ──")

roc_summary <- data.frame(
  Comparison = c("HNSCC vs HC", "HNSCC vs OPMDs", "OPMDs vs HC"),
  n = c(sum(idx_A), sum(idx_B), sum(idx_C)),
  AUC = round(c(as.numeric(ci_A[2]), as.numeric(ci_B[2]), as.numeric(ci_C[2])), 4),
  AUC_CI_lower = round(c(as.numeric(ci_A[1]), as.numeric(ci_B[1]), as.numeric(ci_C[1])), 4),
  AUC_CI_upper = round(c(as.numeric(ci_A[3]), as.numeric(ci_B[3]), as.numeric(ci_C[3])), 4),
  Interpretation = c(
    "Primary diagnostic performance",
    "Disease specificity (HNSCC vs pre-malignant)",
    "Selectivity check (low AUC expected)"
  ),
  stringsAsFactors = FALSE
)
write.csv(roc_summary, "output/tables/Table_cohort3_ROC.csv", row.names = FALSE)
cat("  -> output/tables/Table_cohort3_ROC.csv")

# Prediction scores
pred_scores <- data.frame(
  Sample_ID = cohort3$Sample_ID,
  Group = cohort3$Group,
  Pred_Prob_HNSCC = round(pred_prob, 4)
)
write.csv(pred_scores, "output/tables/Table_cohort3_pred_scores.csv", row.names = FALSE)
cat("  -> output/tables/Table_cohort3_pred_scores.csv")

# ── 9. ROC Figure (3 curves) ─────────────────────────────────────────────

cat("── Step 9: Generating ROC figure ──")

# ROCR objects
rocr_A <- performance(prediction(p_A, y_A), "tpr", "fpr")
rocr_B <- performance(prediction(p_B, y_B), "tpr", "fpr")
rocr_C <- performance(prediction(p_C, y_C), "tpr", "fpr")

COL_A <- "#00A087"  # teal
COL_B <- "#E64B35"  # red
COL_C <- "#4DBBD5"  # blue

pdf("output/figures/Fig_ROC_cohort3.pdf", width = 5, height = 5)
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "False Positive Rate",
     ylab = "True Positive Rate",
     main = "ROC Curves - Cohort 3 Disease Control")
lines(rocr_A@x.values[[1]], rocr_A@y.values[[1]], col = COL_A, lwd = 3)
lines(rocr_B@x.values[[1]], rocr_B@y.values[[1]], col = COL_B, lwd = 3)
lines(rocr_C@x.values[[1]], rocr_C@y.values[[1]], col = COL_C, lwd = 3)
abline(a = 0, b = 1, col = "black", lty = 2, lwd = 1.5)
legend("bottomright",
       legend = c(
         sprintf("HNSCC vs HC (AUC=%.3f)", ci_A[2]),
         sprintf("HNSCC vs OPMDs (AUC=%.3f)", ci_B[2]),
         sprintf("OPMDs vs HC (AUC=%.3f)", ci_C[2])
       ),
       col = c(COL_A, COL_B, COL_C),
       lwd = 3, cex = 0.8, bty = "n")
dev.off()

cat("  -> output/figures/Fig_ROC_cohort3.pdf")

cat("[02_model_validation_cohort3] DONE.")
