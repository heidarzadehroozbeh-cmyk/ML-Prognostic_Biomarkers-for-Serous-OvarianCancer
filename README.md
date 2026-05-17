# ML Prognostic Biomarkers for High-Grade Serous Carcinoma (HGSC)

Machine learning–driven transcriptomic analysis identifying candidate prognostic biomarkers involved in tumor microenvironment (TME) dynamics and progression of **High-Grade Carcionoma (HGSC)**.

---

## Badges

![R](https://img.shields.io/badge/R-%3E%3D4.2-blue)  
![License](https://img.shields.io/badge/license-CC%20BY--NC%204.0-green)  
![Platform](https://img.shields.io/badge/platform-GEO%20GPL570-orange)  
![Status](https://img.shields.io/badge/manuscript-under%20review-yellow)

---

## Overview

This repository contains the computational pipeline used for **integrative transcriptomic analysis and machine learning–based biomarker discovery** in Serous Ovarian Cancer.

The framework integrates multiple independent GEO cohorts and applies a **strict Leave-One-Dataset-Out (LODO) validation strategy** to ensure robust biomarker discovery and cross-cohort generalizability.

The analysis ultimately highlights candidate genes including **SPON1** and **ALDH1A2**, linked to tumor microenvironment remodeling and SOC progression.

---

## Abstract
Serous ovarian cancer (SOC) remains one of the most lethal gynecological malignancies, with limited robust prognostic biomarkers available for clinical use. In this project, we developed a fully reproducible machine learning framework for cross‑cohort transcriptomic analysis of SOC. Publicly available GEO datasets were curated, pre‑processed, and harmonized using standardized normalization and batch correction procedures. Feature selection pipelines combining univariate filtering, penalized regression, and stability selection were applied to identify robust prognostic gene signatures. Multiple machine learning models, including elastic net, random forest, and gradient boosting, were trained and evaluated under rigorous cross‑validation and cross‑cohort validation schemes. Survival modeling and risk stratification were performed to assess the prognostic value of the derived signatures, while pathway‑level analyses and immune deconvolution were used to contextualize the biological relevance of candidate biomarkers. The framework emphasizes methodological transparency, reproducibility, and multi‑cohort robustness, providing a template for biomarker discovery in other cancers. Manuscript submission is in progress; code and workflows will be fully synchronized with the final published version.

## Study Design

The pipeline integrates transcriptomic data across multiple cohorts using the following strategy:

1. GEO dataset acquisition and preprocessing  
2. Differential expression analysis (limma)  
3. Cross-cohort consensus DEG identification  
4. Machine learning model development  
5. Strict cross-study validation (LODO)  
6. Multi-model feature importance ranking  
7. Candidate biomarker selection  
8. Tumor microenvironment (TME) immune module scoring  
9. Linking ML prediction probability with immune system signals

---

## Workflow Diagram

Conceptual pipeline:

```
GEO datasets
     │
     ▼
Data preprocessing
     │
     ▼
Differential expression (limma)
     │
     ▼
Cross‑cohort DEG consensus
     │
     ▼
Machine Learning Models
(ElasticNet / RF / XGBoost)
     │
     ▼
Strict LODO validation
     │
     ▼
Feature ranking
     │
     ▼
Candidate biomarker discovery
     │
     ▼
TME immune module analysis
```

---

## Requirements

Recommended environment:

- **R ≥ 4.2**

Key packages:

```
GEOquery
limma
data.table
metafor
glmnet
randomForest
xgboost
gbm
pROC
tidyverse
survival
ggplot2
```

Install dependencies:

```r
source("scripts/00_setup.R")
```

---

## Quick Start

Example minimal pipeline execution:

```r
# 1. Install dependencies
source("scripts/00_setup.R")

# 2. Download GEO datasets
source("scripts/02_download_geo.R")

# 3. Perform DEG analysis
source("scripts/04_deg_mrna_limma.R")

# 4. Generate consensus genes
source("scripts/05_meta_deg_consensus.R")

# 5. Run ML diagnostic modeling
source("scripts/06_ml_diagnostic_LODO.R")
```

---

## Repository Structure

```
EOC_WNT_TGFb_EMT_Transcriptomics
│
├── data_raw/               # GEO raw files
├── data_processed/         # Processed gene expression matrices
├── meta/                   # Metadata and sample sheets
│
├── results/
│   ├── tables/             # DEG outputs
│   └── plots/              # Standard visualizations
│
├── scripts/
│   ├── 00_setup.R
│   ├── 01_params.R
│   ├── 02_download_geo.R
│   ├── 03_make_sample_sheets.R
│   ├── 04_deg_mrna_limma.R
│   ├── 05_meta_deg_consensus.R
│   ├── 06_ml_diagnostic_LODO.R
│   ├── 07_benchmark_models.R
│   ├── 08_gene_ranking_LODO.R
│   └── 09_select_candidate_biomarkers.R
│
└── README.md
```

---

## Datasets

Platform: **Affymetrix GPL570**

Cohorts used:

- GSE14407  
- GSE38666  
- GSE52037  

DEG filtering criteria:

- Adjusted p-value (FDR) < 0.05  
- |log2FC| > 2  

Consensus genes are defined as **directionally concordant across all cohorts**.

---

## Reproducibility Statement

All analyses were performed using R-based scripts designed for reproducibility across independent datasets. The pipeline emphasizes **cross-cohort validation and strict model generalization** using Leave-One-Dataset-Out training strategies.

---

## Abstract
High-Grade Serous ovarian cancer (HGSC) remains one of the most lethal gynecological malignancies, with limited robust prognostic biomarkers available for clinical application. In this project, we developed a reproducible machine learning framework for cross‑cohort transcriptomic analysis of HGSC.
Publicly available GEO datasets were curated, pre‑processed, and harmonized using standardized normalization and batch‑aware preprocessing procedures. Differential expression analysis was performed using limma, followed by meta‑analytic consensus filtering across datasets.
Feature selection pipelines combining univariate filtering, penalized regression, and stability‑based ranking were applied to identify robust candidate genes. Multiple machine learning models—including Elastic Net, Random Forest, Gradient Boosting, and XGBoost—were trained and evaluated using a strict Leave‑One‑Dataset‑Out (LODO) validation framework.
Model performance was evaluated using ROC‑based metrics and cross‑cohort predictive consistency. Candidate biomarkers were further contextualized through immune module scoring to investigate associations between ML prediction probability and tumor microenvironment activity.
This framework emphasizes methodological transparency, reproducibility, and cross‑cohort robustness, providing a generalizable strategy for transcriptomic biomarker discovery in cancer research.n.

---

## Author

Dr Roozbeh Heidarzadehpilehrood  
Human Genetics, Genomics & Transcriptomics

Email  
roozbeh.heidarzadeh@gmail.com  
heidarzadeh.roozbeh@gmail.com

---

## 📄 Citation

If you use this pipeline, please cite:

Heidarzadehpilehrood R, et.al (2026)
.......................
Repository DOI  
.......................

---

## 📝 License

This project is licensed under:

**Creative Commons CC BY‑NC 4.0**

You may share and adapt the material for **non‑commercial use** with proper attribution.


---
