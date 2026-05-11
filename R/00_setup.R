#=============================================================================
# 00_setup.R — Environment and dependencies
#=============================================================================
# Usage: source("R/00_setup.R") at the top of each analysis script
#=============================================================================

# Set working directory to project root (if running interactively)
if (interactive()) {
  # Uncomment and modify if needed:
  # setwd("/path/to/HNSCC-Metabolomics")
}

# Required packages
required_packages <- c(
  "randomForest",  # Random forest classifier
  "pROC",          # ROC analysis and CI
  "ROCR",          # ROC curve plotting
  "ggplot2"        # Figures
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Create output directories
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

cat("Environment ready.")
cat(sprintf("R version: %s", R.version.string))
cat(sprintf("randomForest: %s", packageVersion("randomForest")))
cat(sprintf("pROC: %s", packageVersion("pROC")))
