# ML-Prognostic_Biomarkers-for-Serous-OvarianCancer
Developed ML algorisms on transcriptomics landscapes identi-fies SPON1 and ALDH1A2 as candidate prognostic biomarkers in TME and serous ovarian cancer progression

---

## Table of Contents

1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Repository Structure](#repository-structure)
4. [Usage](#usage)
5. [Scripts Description](#scripts-description)
6. [Datasets](#datasets)
7. [Author](#author)
8. [Citation](#citation)
9. [License](#license)

---

## Overview

This project integrates multiple transcriptomic datasets of SOC to identify Candidate gene panel and validate key and core genes. The workflow includes:

* Data mining from GEO datasets
* Expression preprocessing
* Cohort-wise differential expression
* Pathway analysis
* Cross-cohort meta-integration
* Strict leave-one-dataset-out (LODO) diagnostic machine learning
* Multi-model feature ranking
* Validation of Candidate Genes
* Tumor microenvironment (TME) assessment
* Linking ML predicted probability to immune modules

---

## Requirements

*   **R:** Version ≥ 4.2 recommended.
*   **R Packages:**
    *   `GEOquery`
    *   `limma`
    *   `data.table`
    *   `metafor`
    *   `glmnet`
    *   `pROC`
    *   `randomForest`
    *   `gbm`
    *   `xgboost`
    *   `tidyverse`
    *   `survival` (for survival analysis)
    *   `ggplot2` (for plotting)
    *   Bioconductor annotation packages (e.g., `hgu133plus2.db`)

All packages can be installed automatically using `scripts/00_setup.R`.

---

## Repository Structure

```
EOC_WNT_TGFb_EMT_Transcriptomics/
│
├── data_raw/                                               # Raw GEO files (downloaded GSEMatrix objects)
├── data_processed/                                         # Processed gene-level expression matrices
├── meta/                                                   # Sample sheets and metadata CSVs
├── results/
│   ├── tables/                                             # DEG results per dataset
│   └── plots/                                              # Volcano plots per dataset
├── scripts/
│   ├── 00_setup.R                                          # Install required packages
│   ├── 01_params.R                                         # Project parameters, DEG thresholds, paths
│   ├── 02_download_geo.R                                   # Download GEO GSEMatrix objects
│   ├── 03_make_sample_sheets.R                             # Create/validate sample sheets
│   ├── 04_deg_mrna_limma.R                                 # Differential expression analysis (limma)
│   ├── 05_meta_deg_consensus.R                             # Create meta_deg_consensus
│   ├── 06_ml_diagnostic_LODO.R                             # Diagnostic ML (EOC vs Normal) with cross-study LODO validation
│   ├── 07_DEGs_expression_Panel.R                          # Figure Panels
│   ├── 08_Predictions_confusion.R                          # Create meta_deg_consensus
│   ├── 09_figures_STRICT_ROC_and_expression_panels.R       # Publication-grade ROC + expression panels
│   ├── 10_benchmark_models_STRICT_LODO.R                   # Models: ElasticNet (glmnet), RandomForest (ranger), XGBoost (xgboost)
│   ├── 11_multimodel_gene_ranking_STRICT_LODO_top500.R     # Multi-model benchmarking + gene ranking on TRAIN significant genes within each LODO fold
│   ├── 12_select_candidate_biomarkers_STRICT_LODO.R        # Candidate biomarker selection using CLASSIC ML models (ElasticNet / SVM / RF / XGBoost / GBM)
│   ├── 13_BIGPANEL_STRICT_candidates_parts.R               # BigPanel Outputs
│   ├── 14_FINAL22_multigene_expression_Q1.R                # 22 FINAL genes (Up/Down) with Box Violin Points
│   ├── 15_cripts/16_ALGO_training_curves_COREset22.R

           # BigPanel Outputs
│   ├── 16_BIGPANEL_STRICT_candidates_parts.R               # BigPanel Outputs
│   ├── 16_BIGPANEL_STRICT_candidates_parts.R               # BigPanel Outputs
│   ├── 16_BIGPANEL_STRICT_candidates_parts.R               # BigPanel Outputs


# Inputs (from script 11):
└── README.md                 # Project description (this file)
```

---

## Datasets

* **mRNA datasets (GPL570)**: `GSE14407`, `GSE38666`, `GSE52037`

DEG thresholds used:

* **Adjusted P-value (FDR) < 0.05**
* **|log2 Fold Change| > 2** (for mRNA)

Consensus DEGs are defined as **directionally concordant across all 3 mRNA datasets**.

---

## Author

Dr Roozbeh Heidarzadehpilehrood
Affiliation: Independent researcher, Human Genetics, Genomics & Transcriptomics
Contact: heidarzadeh.roozbeh@gmail.com

---

## Citation

If you use this code or parts of the pipeline, please cite both:

> Heidarzadehpilehrood R, Ling K-H, Abdul Hamid H (2026) Integrative transcriptomic analysis of WNT/TGFβ-driven EMT pathways and drug-gene interaction networks in epithelial ovarian cancer. Advances in Cancer Biology - Metastasis 16:100178. [https://doi.org/10.1016/j.adcanc.2026.100178].

> Heidarzadehpilehrood R, 2026. GitHub repository: EOC_WNT_TGFb_EMT_Transcriptomics_2026. (https://zenodo.org/records/18711967)


---

## License

This project is released under **CC BY-NC 4.0**.

You are free to **share** and **adapt** the material for **non-commercial purposes**, provided appropriate credit is given and modifications are indicated.

---
