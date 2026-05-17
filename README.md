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
