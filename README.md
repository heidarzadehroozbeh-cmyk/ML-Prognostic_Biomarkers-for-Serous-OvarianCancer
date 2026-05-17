# ML Prognostic Biomarkers for Serous Ovarian Cancer

Machine learning–driven transcriptomic analysis identifying candidate prognostic biomarkers involved in tumor microenvironment (TME) dynamics and progression of **Serous Ovarian Cancer (SOC)**.

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
│   ├── 10_benchmark_models.R
│   ├── 11_gene_ranking_LODO.R
│   ├── 12_select_candidate_biomarkers.R
│   ├── 21_immune_scoring.R
│   └── 22_ML_vs_Immune.R
│
└── README.md
```

Note: **High‑resolution figure generation scripts for the manuscript are intentionally excluded and will be released after publication.**

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

## Manuscript Status

This repository accompanies a research manuscript currently **under peer review**.  
Some scripts related to **publication-grade figure generation** are temporarily withheld and will be released after formal publication.

---

## Author

Dr Roozbeh Heidarzadehpilehrood  
Independent Researcher — Human Genetics, Genomics & Transcriptomics

Email  
roozbeh.heidarzadeh@gmail.com  
heidarzadeh.roozbeh@gmail.com

---

## Citation

If you use this pipeline, please cite:

Heidarzadehpilehrood R, et.al (2026)
.......................
Repository DOI  
.......................

---

## License

This project is licensed under:

**Creative Commons CC BY‑NC 4.0**

You may share and adapt the material for **non‑commercial use** with proper attribution.






```markdown
# ML Prognostic Biomarkers for Serous Ovarian Cancer

**Machine learning–driven transcriptomic analysis for identifying robust prognostic biomarkers associated with tumor microenvironment (TME) dynamics and progression of Serous Ovarian Cancer (SOC).**

---

## 📚 Overview

This repository contains a reproducible computational pipeline for **integrative transcriptomic analysis and machine learning–based biomarker discovery** in serous ovarian cancer. The framework utilizes multiple GEO datasets and implements a strict **Leave-One-Dataset-Out (LODO)** validation strategy.

---

## ⚙ Repository Structure

```plaintext
ML-Prognostic_Biomarkers-for-Serous-OvarianCancer
│
├── scripts/                            # Analysis scripts
│   ├── 00_setup.R                      # Install dependencies
│   ├── 01_params.R                     # Global parameters
│   ├── 02_download_geo.R               # Download GEO datasets
│   ├── 03_make_sample_sheets.R         # Construct sample metadata
│   ├── 04_deg_mrna_limma.R             # DEG analysis using limma
│   ├── 05_meta_deg_consensus.R         # Cross-cohort consensus DEGs
│   ├── 06_ml_diagnostic_LODO.R         # ML training with LODO
│   ├── 07_benchmark_models_STRICT_LODO.R# Benchmarking ML models
│   ├── 08_multimodel_gene_ranking_STRICT_LODO_top500.R # Feature ranking
│   └── 09_select_candidate_biomarkers_STRICT_LODO.R  # Candidate selection
│
├── README.md                           # Project documentation
└── LICENSE                             # Licensing information
```

---

## 📊 Study Design

1. **Acquire GEO datasets**  
2. **Normalize and preprocess** datasets  
3. **Perform differential expression** analysis using limma  
4. **Generate cross-cohort consensus DEGs**  
5. **Train machine learning models** with strict LODO validation  
6. **Benchmark multiple algorithms**  
7. **Rank genes** based on feature importance  
8. **Select robust biomarkers**  
9. **Evaluate associations** with immune signals  

---

## 🚀 Quick Start

### Step-by-step Instructions

```r
# Install dependencies
source("scripts/00_setup.R")

# Download and preprocess GEO data
source("scripts/02_download_geo.R")

# Construct sample metadata
source("scripts/03_make_sample_sheets.R")

# Execute DEG analysis
source("scripts/04_deg_mrna_limma.R")

# Generate consensus DEGs
source("scripts/05_meta_deg_consensus.R")

# Train ML models with LODO validation
source("scripts/06_ml_diagnostic_LODO.R")

# Benchmark models
source("scripts/07_benchmark_models_STRICT_LODO.R")

# Rank top genes
source("scripts/08_multimodel_gene_ranking_STRICT_LODO_top500.R")

# Select candidate biomarkers
source("scripts/09_select_candidate_biomarkers_STRICT_LODO.R")
```

---

## 📈 Example Results

![Example Workflow](path_to_your_workflow_image)  
*Sample visualization of the workflow and results.*

---

## 🔗 Additional Information

- Datasets analyzed include: **GSE14407**, **GSE38666**, **GSE52037**.
- Repository is currently **under peer review** and will be updated post-publication.
- For any inquiries or collaboration opportunities, contact:  
  Dr. Roozbeh Heidarzadehpilehrood  
  [Email](mailto:roozbeh.heidarzadeh@gmail.com)

---

## 📝 License

Creative Commons **CC BY-NC 4.0** — non-commercial use allowed with attribution.

---

## 📄 Citation

If you use this repository, please cite:  
Heidarzadehpilehrood R. et al. (2026). Machine learning–based discovery of prognostic biomarkers in serous ovarian cancer.

Repository DOI: *To be assigned via Zenodo*
```
