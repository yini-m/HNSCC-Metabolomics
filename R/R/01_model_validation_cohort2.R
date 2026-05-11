#=============================================================================
# 01_model_validation_cohort2.R
# Cohort 2 External Validation — 8-Metabolite RF Model
#=============================================================================
# Author:   Hexin Ma, Xinran Zhao
# Date:     2026-05
#
# Purpose:  Validates the 8-metabolite random forest model on the independent
#           external cohort (Cohort 2, n=171) using per-cohort independent
#           Z-score normalisation. Reports overall AUC with 95% CI.
#
# Input:    data/cohort1_8metabolites.csv
#           data/cohort2_8metabolites.csv
#           data/selected_features.csv
#
# Output:   output/tables/Table_cohort2_ROC.csv
#           output/tables/Table_cohort2_pred_scores.csv
#           output/figures/Fig_ROC_cohort2.pdf
#
# Usage:    Rscript R/01_model_validation_cohort2.R   (from project root)
# License:  GPL-3.0
#=============================================================================

source("R/00_setup.R")

library(randomForest)
library(pROC)
library(ROCR)

cat("=============================================================")
cat("01  Cohort 2 External Validation")
cat("=============================================================")

# ── 1. Load Data ─────────────────────────────────────────────────────────

cat("── Step 1: Loading data ──")

features <- read.csv("data/selected_features.csv")
metabolites <- features$Feature

cohort1 <- read.csv("data/cohort1_8metabolites.csv")
cohort2 <- read.csv("data/cohort2_8metabolites.csv")

cat(sprintf("  Cohort 1 (training): %d samples", nrow(cohort1)))
cat(sprintf("  Cohort 2 (external): %d samples (T=%d, N=%d)",
            nrow(cohort2),
            sum(cohort2$Group == "T"),
            sum(cohort2$Group == "N")))

# ── 2. Independent Z-Score Normalisation ──────────────────────────────────

cat("── Step 2: Independent Z-score normalisation ──")

# Cohort 1: Z-score using its own parameters
for (m in metabolites) {
  mu <- mean(cohort1[[m]], na.rm = TRUE)
  s  <- sd(cohort1[[m]], na.rm = TRUE)
  cohort1[[m]] <- (cohort1[[m]] - mu) / s
}

# Cohort 2: Z-score using its own parameters
for (m in metabolites) {
  mu <- mean(cohort2[[m]], na.rm = TRUE)
  s  <- sd(cohort2[[m]], na.rm = TRUE)
  cohort2[[m]] <- (cohort2[[m]] - mu) / s
}

cat("  Cohort 1 post-Z mean (should be ~0):")
print(round(colMeans(cohort1[, metabolites]), 3))
cat("  Cohort 2 post-Z mean (should be ~0):")
print(round(colMeans(cohort2[, metabolites]), 3))

# ── 3. Train RF on Cohort 1 ──────────────────────────────────────────────

cat("── Step 3: Training RF on Cohort 1 (Z-scored) ──")

y_train <- factor(ifelse(cohort1$Group == "T", "HNSCC", "HC"),
                  levels = c("HC", "HNSCC"))

train_df <- cohort1[, metabolites, drop = FALSE]
train_df$Label <- y_train

set.seed(42)
rf_model <- randomForest(Label ~ ., data = train_df, ntree = 500, importance = TRUE)
cat("  RF trained (ntree=500)")

# ── 4. Predict Cohort 2 ──────────────────────────────────────────────────

cat("── Step 4: Predicting Cohort 2 ──")

test_df <- cohort2[, metabolites, drop = FALSE]
pred_prob <- predict(rf_model, newdata = test_df, type = "prob")[, "HNSCC"]

y_true <- ifelse(cohort2$Group == "T", 1, 0)

# ── 5. ROC Analysis ──────────────────────────────────────────────────────

cat("── Step 5: ROC analysis ──")

roc_c2 <- roc(y_true, pred_prob, levels = c(0, 1), direction = "<", quiet = TRUE)
ci_c2  <- ci.auc(roc_c2)
best_c2 <- coords(roc_c2, x = "best", best.method = "youden",
                  ret = c("threshold", "sensitivity", "specificity",
                          "ppv", "npv", "accuracy"))

cat(sprintf("  AUC = %.4f (95%% CI: %.4f-%.4f)",
            ci_c2[2], ci_c2[1], ci_c2[3]))
cat(sprintf("  Sensitivity = %.1f%%", best_c2$sensitivity * 100))
cat(sprintf("  Specificity = %.1f%%", best_c2$specificity * 100))
cat(sprintf("  Accuracy = %.1f%%", best_c2$accuracy * 100))

# Export ROC summary
roc_df <- data.frame(
  Cohort = "Cohort 2",
  n = nrow(cohort2),
  n_HNSCC = sum(y_true == 1),
  n_HC = sum(y_true == 0),
  AUC = round(as.numeric(ci_c2[2]), 4),
  AUC_CI_lower = round(as.numeric(ci_c2[1]), 4),
  AUC_CI_upper = round(as.numeric(ci_c2[3]), 4),
  Sensitivity = round(best_c2$sensitivity * 100, 1),
  Specificity = round(best_c2$specificity * 100, 1),
  PPV = round(best_c2$ppv * 100, 1),
  NPV = round(best_c2$npv * 100, 1),
  Accuracy = round(best_c2$accuracy * 100, 1)
)
write.csv(roc_df, "output/tables/Table_cohort2_ROC.csv", row.names = FALSE)
cat("  -> output/tables/Table_cohort2_ROC.csv")

# Export prediction scores
pred_scores <- data.frame(
  Sample_ID = cohort2$Sample_ID,
  Group = cohort2$Group,
  Pred_Prob_HNSCC = round(pred_prob, 4)
)
write.csv(pred_scores, "output/tables/Table_cohort2_pred_scores.csv", row.names = FALSE)
cat("  -> output/tables/Table_cohort2_pred_scores.csv")

# ── 6. ROC Figure ─────────────────────────────────────────────────────────

cat("── Step 6: Generating ROC figure ──")

rocr_pred <- prediction(pred_prob, y_true)
rocr_perf <- performance(rocr_pred, "tpr", "fpr")

pdf("output/figures/Fig_ROC_cohort2.pdf", width = 5, height = 5)
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "False Positive Rate",
     ylab = "True Positive Rate",
     main = "ROC Curve - Cohort 2 External Validation")
lines(rocr_perf@x.values[[1]], rocr_perf@y.values[[1]],
      col = "#00A087", lwd = 3)
abline(a = 0, b = 1, col = "black", lty = 2, lwd = 1.5)
legend("bottomright",
       legend = sprintf("AUC = %.3f (95%% CI: %.3f-%.3f)",
                        ci_c2[2], ci_c2[1], ci_c2[3]),
       col = "#00A087", lwd = 3, cex = 0.85, bty = "n")
dev.off()

cat("  -> output/figures/Fig_ROC_cohort2.pdf")

cat("[01_model_validation_cohort2] DONE.")
