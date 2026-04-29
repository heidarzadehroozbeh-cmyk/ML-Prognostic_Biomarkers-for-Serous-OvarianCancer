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
├── data_raw/                 # Raw GEO files (downloaded GSEMatrix objects)
├── data_processed/           # Processed gene-level expression matrices
├── meta/                     # Sample sheets and metadata CSVs
├── results/
│   ├── tables/               # DEG results per dataset
│   └── plots/                # Volcano plots per dataset
├── scripts/
│   ├── 00_setup.R            # Install required packages
│   ├── 01_params.R           # Project parameters, DEG thresholds, paths
│   ├── 02_download_geo.R     # Download GEO GSEMatrix objects
│   ├── 03_make_sample_sheets.R # Create/validate sample sheets
│   ├── 04_deg_mrna_limma.R   # Differential expression analysis (limma)
│   └── 05_Volcano Plots.R    # Generate volcano plots for DEGs
└── README.md                 # Project description (this file)
```

---

## Usage

1. **Setup environment**

   ```R
   source("scripts/00_setup.R")
   ```

2. **Set project parameters** (optional modifications)

   ```R
   source("scripts/01_params.R")
   ```

3. **Download GEO datasets**

   ```R
   source("scripts/02_download_geo.R")
   ```

4. **Create and validate sample sheets**

   ```R
   source("scripts/03_make_sample_sheets.R")
   ```

5. **Run differential expression analysis (mRNA)**

   ```R
   source("scripts/04_deg_mrna_limma.R")
   ```

6. **Generate volcano plots**

   ```R
   source("scripts/05_Volcano Plots.R")
   ```

---

## Scripts Description

| Script                    | Description                                                                                                       |
| ------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| `00_setup.R`              | Installs required CRAN and Bioconductor packages.                                                                 |
| `01_params.R`             | Defines datasets, expected sample counts, DEG thresholds, and project folder paths.                               |
| `02_download_geo.R`       | Downloads GEO GSEMatrix objects and saves metadata for each dataset.                                              |
| `03_make_sample_sheets.R` | Creates or validates sample sheets for each dataset; ensures `EOC` and `Normal` counts match expectations.        |
| `04_deg_mrna_limma.R`     | Performs differential expression analysis using **limma**; collapses probes to gene symbols and saves DEG tables. |
| `05_Volcano Plots.R`      | Generates volcano plots for each dataset, highlighting DEGs above thresholds and labeling top genes.              |

---

## Datasets

* **mRNA datasets (GPL570)**: `GSE14407`, `GSE38666`, `GSE52037`
* **miRNA dataset (GPL20712)**: `GSE216150`

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
