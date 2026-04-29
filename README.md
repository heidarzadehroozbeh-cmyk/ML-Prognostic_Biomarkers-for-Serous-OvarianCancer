# ML-Prognostic-Biomarkers-SOC
### Integrative Machine Learning Pipeline for Transcriptomic Landscape Analysis in Serous Ovarian Cancer

This repository contains the computational pipeline used to identify and validate candidate prognostic biomarkers in the Tumor Microenvironment (TME) of Serous Ovarian Cancer (SOC). The workflow focuses on multi-cohort integration and robust Machine Learning validation.

---

## 📋 Table of Contents
- [Overview](#overview)
- [Methodological Workflow](#methodological-workflow)
- [Requirements](#requirements)
- [Repository Structure](#repository-structure)
- [Datasets](#datasets)
- [Citation](#citation)
- [License](#license)

---

## 🔍 Overview
This project integrates multiple SOC transcriptomic datasets to identify a robust candidate gene panel. The core of the study is a **Strict Leave-One-Dataset-Out (LODO)** diagnostic framework, ensuring that the identified biomarkers (including **SPON1** and **ALDH1A2**) maintain high predictive power across independent cohorts and diverse clinical settings.

---

## ⚙️ Methodological Workflow
1.  **Data Harmonization:** Mining and preprocessing GEO datasets (GPL570).
2.  **Meta-Analysis:** Cohort-wise differential expression and cross-cohort consensus integration.
3.  **Machine Learning:** Diagnostic modeling using ElasticNet, Random Forest, and XGBoost.
4.  **Robust Validation:** Multi-model feature ranking and strict LODO cross-study validation.
5.  **Systems Biology:** Linking ML-predicted probabilities to immune modules and TME assessment.

---

## 💻 Requirements
- **R Version:** ≥ 4.2
- **Key Packages:** `limma`, `metafor`, `glmnet`, `randomForest`, `xgboost`, `tidyverse`, `survival`, `pROC`.

You can initialize the environment by running:
```R
source("scripts/00_setup.R")
```

---

## 📁 Repository Structure (Core Pipeline)
*Note: High-resolution visualization scripts for the final manuscript figures are currently restricted and will be fully released upon publication.*

```text
scripts/
├── 00_setup.R                    # Environment setup & dependencies
├── 01_params.R                   # Global parameters & DEG thresholds
├── 02_download_geo.R             # Automated data acquisition from GEO
├── 03_make_sample_sheets.R       # Metadata validation
├── 04_deg_mrna_limma.R           # Differential Expression Analysis
├── 05_meta_deg_consensus.R       # Consensus Meta-integration
├── 06_ml_diagnostic_LODO.R       # Primary ML modeling (LODO framework)
├── 10_benchmark_models.R         # Model comparison (RF, XGBoost, ENET)
├── 11_gene_ranking_LODO.R        # Multi-model feature importance ranking
├── 12_select_biomarkers.R        # Identification of candidate biomarkers
├── 21_immune_scoring.R           # TME immune module scoring (ssGSEA-lite)
├── 22_ML_vs_Immune.R             # Correlation of ML outputs & Biology
└── 23_volcano_standard.R         # Standard DEG visualization
```

---

## 📊 Datasets
The analysis utilizes mRNA expression profiles from the **GPL570** platform:
- **GSE14407**, **GSE38666**, **GSE52037**
- **Inclusion Criteria:** Adjusted P-value (FDR) < 0.05 & |log₂FC| > 2.

---

## ✍️ Author
**Dr. Roozbeh Heidarzadeh-Pilehrood**  
*Independent Researcher | Human Genetics & Transcriptomics*  
📧 [roozbeh.heidarzadeh@gmail.com](mailto:roozbeh.heidarzadeh@gmail.com)

---

## 📜 Citation
If you use this pipeline, please cite:

1.  **Heidarzadehpilehrood R**, et al. (2026). *Integrative transcriptomic analysis of WNT/TGFβ-driven EMT pathways...* **Advances in Cancer Biology - Metastasis**, 16:100178. [DOI: 10.1016/j.adcanc.2026.100178]
2.  **Heidarzadehpilehrood R**, (2026). GitHub: *EOC_WNT_TGFb_EMT_Transcriptomics_2026*. [Zenodo](https://zenodo.org/records/18711967)

---

## 📄 License
This project is licensed under **CC BY-NC 4.0**. (Non-commercial use with attribution).

---
