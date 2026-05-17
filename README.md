# ML Prognostic Biomarkers for High-Grade Serous Carcinoma (HGSC)

Machine learning–driven transcriptomic analysis identifying candidate prognostic biomarkers involved in tumor microenvironment (TME) dynamics and progression of **High-Grade Serous Carcionoma (HGSC)**.

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
High-Grade Serous Carcionoma (HGSC) remains one of the most lethal gynecological malignancies, with limited robust prognostic biomarkers available for clinical use. In this project, we developed a fully reproducible machine learning framework for cross‑cohort transcriptomic analysis of HGSC. Publicly available GEO datasets were curated, pre‑processed, and harmonized using standardized normalization and batch correction procedures. Feature selection pipelines combining univariate filtering, penalized regression, and stability selection were applied to identify robust prognostic gene signatures. Multiple machine learning models, including elastic net, random forest, and gradient boosting, were trained and evaluated under rigorous cross‑validation and cross‑cohort validation schemes. Survival modeling and risk stratification were performed to assess the prognostic value of the derived signatures, while pathway‑level analyses and immune deconvolution were used to contextualize the biological relevance of candidate biomarkers. The framework emphasizes methodological transparency, reproducibility, and multi‑cohort robustness, providing a template for biomarker discovery in other cancers. Manuscript submission is in progress; code and workflows will be fully synchronized with the final published version.

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


This project is licensed under the **MIT License**.

Copyright (c) 2026 Roozbeh Heidarzadehpilehrood






# 🧬 Machine Learning Prognostic Biomarkers for High‑Grade Serous Carcinoma (HGSC)

Machine learning–driven transcriptomic biomarker discovery framework for identifying prognostic genes involved in tumor microenvironment (TME) remodeling and disease progression in High‑Grade Serous Ovarian Carcinoma (HGSC).

🧬 Human Genomics • 🧫 Transcriptomics • 🧪 Tumor Microenvironment • 🤖 Machine Learning • 📊 Survival Modeling

---

## 🏷 Badges

![R](https://img.shields.io/badge/R-%3E%3D4.2-blue)
![Machine Learning](https://img.shields.io/badge/Machine%20Learning-ElasticNet%20%7C%20RF%20%7C%20XGBoost-purple)
![Platform](https://img.shields.io/badge/platform-GEO%20GPL570-orange)
![Transcriptomics](https://img.shields.io/badge/data-Transcriptomics-blue)
![Reproducibility](https://img.shields.io/badge/Reproducible-Yes-brightgreen)
![Status](https://img.shields.io/badge/manuscript-under%20review-yellow)
![License](https://img.shields.io/badge/license-MIT-green)

---

# 🔬 Overview

This repository contains a fully reproducible computational framework for cross‑cohort transcriptomic biomarker discovery in High‑Grade Serous Ovarian Carcinoma (HGSC).

The pipeline integrates multiple independent GEO transcriptomic cohorts and applies a strict Leave‑One‑Dataset‑Out (LODO) validation strategy to ensure robust feature selection and cross‑study generalization.

The analytical workflow combines:

• 🧬 Differential expression meta‑analysis  
• 🤖 Machine learning‑based feature selection  
• 📊 Cross‑cohort model validation  
• 🧪 Tumor microenvironment immune analysis  
• 📈 Survival‑linked biomarker prioritization  

Candidate biomarkers emerging from this framework include genes such as:

SPON1  
ALDH1A2

These genes are associated with tumor microenvironment remodeling and ovarian cancer progression.

---

# 📑 Abstract

High‑Grade Serous Ovarian Carcinoma (HGSC) remains one of the most lethal gynecologic malignancies, largely due to late diagnosis and the absence of robust prognostic biomarkers suitable for clinical translation.

In this study we developed a cross‑cohort transcriptomic biomarker discovery framework integrating multiple independent GEO datasets. Raw microarray datasets were systematically curated, preprocessed, normalized, and harmonized using standardized pipelines to minimize technical heterogeneity.

Differential gene expression analysis was performed using the limma empirical Bayes framework, followed by cross‑study consensus DEG identification across independent cohorts.

To identify robust predictive biomarkers we implemented a machine learning ensemble framework combining:

Elastic Net regularized regression  
Random Forest  
Gradient Boosting Machines  
XGBoost  

Model robustness was evaluated using a strict Leave‑One‑Dataset‑Out validation strategy, ensuring that candidate biomarkers demonstrate consistent predictive performance across independent cohorts.

To contextualize biological significance additional analyses were performed including:

Tumor microenvironment immune module scoring  
Pathway enrichment analysis  
Survival‑based risk stratification  

The pipeline emphasizes methodological transparency, reproducibility, and cross‑cohort robustness.

Manuscript submission is currently in progress.

---

# 🧠 Study Design

The integrative analysis pipeline follows the workflow below:

1. 🧬 GEO dataset acquisition and curation  
2. 🧹 Data preprocessing and normalization  
3. 📊 Differential expression analysis (limma)  
4. 🔗 Cross‑cohort consensus DEG identification  
5. 🤖 Machine learning model training  
6. 🔁 Strict Leave‑One‑Dataset‑Out validation  
7. 📈 Model benchmarking and performance comparison  
8. 🧬 Feature importance ranking  
9. 🧪 Tumor microenvironment immune module scoring  
10. 🎯 Candidate biomarker prioritization  

---

# 🔁 Workflow Diagram

GEO datasets  
     │  
     ▼  
Data preprocessing  
     │  
     ▼  
Differential expression analysis (limma)  
     │  
     ▼  
Cross‑cohort consensus DEG identification  
     │  
     ▼  
Machine Learning Models  
ElasticNet | Random Forest | XGBoost | GBM  
     │  
     ▼  
Strict LODO validation  
     │  
     ▼  
Feature importance ranking  
     │  
     ▼  
Candidate biomarker discovery  
     │  
     ▼  
Tumor Microenvironment analysis  
     │  
     ▼  
Survival‑linked biomarker prioritization  

---

# ⚙️ Requirements

Recommended environment:

R ≥ 4.2

---

## 📦 Core Bioinformatics Packages

GEOquery  
limma  
edgeR  
Biobase  
BiocGenerics  
SummarizedExperiment  
affy  
oligo  
annotate  
hgu133plus2.db  
AnnotationDbi  

---

## 📊 Data Processing

data.table  
dplyr  
tidyr  
tibble  
readr  
stringr  
forcats  
purrr  
janitor  

---

## 🤖 Machine Learning

glmnet  
caret  
randomForest  
xgboost  
gbm  
e1071  
ranger  
nnet  
kernlab  
MLmetrics  

---

## 📈 Statistics

survival  
survminer  
metafor  
pROC  
ROCR  
boot  

---

## 🧬 Functional Analysis

clusterProfiler  
fgsea  
ReactomePA  
DOSE  
enrichR  
GSVA  

---

## 🧪 Tumor Microenvironment

GSVA  
estimate  
immunedeconv  
MCPcounter  
xCell  
CIBERSORT  

---

## 📊 Visualization

ggplot2  
ComplexHeatmap  
pheatmap  
corrplot  
cowplot  
patchwork  
RColorBrewer  
viridis  

---

## 🧰 Utility

optparse  
argparse  
here  
glue  
fs  
yaml  

---

## Install dependencies

source("scripts/00_setup.R")

---

# 🚀 Quick Start

Example minimal pipeline execution:

source("scripts/00_setup.R")

source("scripts/02_download_geo.R")

source("scripts/04_deg_mrna_limma.R")

source("scripts/05_meta_deg_consensus.R")

source("scripts/06_ml_diagnostic_LODO.R")

---

# 📂 Repository Structure

ML-Prognostic_Biomarkers-for-Serous-OvarianCancer

data_raw  
data_processed  
meta  

scripts  
00_setup.R  
01_params.R  
02_download_geo.R  
03_make_sample_sheets.R  
04_deg_mrna_limma.R  
05_meta_deg_consensus.R  
06_ml_diagnostic_LODO.R  
07_benchmark_models.R  
08_gene_ranking_LODO.R  
09_select_candidate_biomarkers.R  

README.md

---

# 🧬 Datasets

Platform:

Affymetrix Human Genome U133 Plus 2.0 (GPL570)

GEO cohorts:

GSE14407  
GSE38666  
GSE52037  

Filtering criteria:

FDR < 0.05  
|log2FC| > 2  

Consensus genes are defined as directionally concordant across all cohorts.

---

# 🔁 Reproducibility

All analyses are implemented using fully scripted R workflows designed for reproducible cross‑cohort analysis.

The pipeline emphasizes:

reproducible preprocessing  
transparent statistical testing  
cross‑dataset validation  
machine learning robustness  

The Leave‑One‑Dataset‑Out validation strategy ensures that biomarkers generalize across independent cohorts.

---

# 👨‍🔬 Author

Dr Roozbeh Heidarzadehpilehrood

Human Genetics • Genomics • Transcriptomics • Computational Biology

Email

roozbeh.heidarzadeh@gmail.com  
heidarzadeh.roozbeh@gmail.com  

---

# 📄 Citation

If you use this pipeline please cite:

Heidarzadehpilehrood R et al. (2026)

Machine learning‑driven transcriptomic biomarker discovery in HGSC.

Repository DOI will be updated after Zenodo release.

---

# 📜 License

MIT License

Copyright (c) 2026  
Roozbeh Heidarzadehpilehrood

See LICENSE file for details.


See the [LICENSE](LICENSE) file for the full license text.


---
