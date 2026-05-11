# HNSCC-Metabolomics
Serum metabolomics-based diagnostic model for head and neck squamous cell carcinoma (HNSCC)

## Overview

This repository contains the core analysis code for validating an 8-metabolite random forest diagnostic model on two independent cohorts, as described in:


## Repository Structure

HNSCC-Metabolomics/
├── README.md
├── LICENSE
├── R/
│   ├── 00_setup.R                          # Environment and dependencies
│   ├── 01_model_validation_cohort2.R       # External validation (Cohort 2)
│   └── 02_model_validation_cohort3.R       # Disease control validation (Cohort 3)
├── data/
│   ├── cohort1_8metabolites.csv            # Discovery cohort (n=498, training)
│   ├── cohort2_8metabolites.csv            # External validation (n=171)
│   ├── cohort3_8metabolites.csv            # Disease control (n=269)
│   └── selected_features.csv              # 8 diagnostic metabolites
└── output/
    ├── tables/                             # Generated CSV results
    └── figures/                            # Generated PDF figures

## Requirements

- R >= 4.2.0
- randomForest >= 4.7
- pROC >= 1.18
- ROCR >= 1.0
- ggplot2 >= 3.4

## How to Use

From project root:
Rscript R/01_model_validation_cohort2.R
Rscript R/02_model_validation_cohort3.R

Each script independently loads Cohort 1 data, trains the RF model, and validates on the target cohort.

## Data Description

All data files contain scaled concentrations of 8 diagnostic metabolites. Per-cohort independent Z-score normalisation is performed within each analysis script. Sample IDs are anonymised.


## Citation

If you use this code or data, please cite:

> Ma H, Zhao X, et al. Serum metabolomics-based diagnostic model for head and neck squamous cell carcinoma: a multicenter validation study. *eBioMedicine*. 2026. (Under Review)

## License

GPL-3.0
