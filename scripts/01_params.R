# scripts/01_params.R
set.seed(123)

# --- Datasets from manuscript ---
GSE_mrna <- c("GSE14407", "GSE38666", "GSE52037")  # GPL570

expected_counts <- list(
  GSE14407 = c(HGSC = 12, Normal = 12),
  GSE38666 = c(HGSC = 13, Normal = 12),
  GSE52037 = c(HGSC = 9, Normal = 10)
)

# --- DEG thresholds (as reported) ---
DEG_FDR <- 0.05
DEG_abs_logFC <- 2   # |log2FC| > 2 used in your results/figs

# --- Consensus definition ---
CONSENSUS_directional <- TRUE  # directionally concordant across all 3 mRNA sets

# --- Paths ---
DIR_RAW  <- "data_raw"
DIR_GEO  <- file.path(DIR_RAW, "GEO")
DIR_PROC <- "data_processed"
DIR_META <- "meta"
DIR_RES  <- "results"
DIR_TAB  <- file.path(DIR_RES, "tables")
DIR_MOD  <- file.path(DIR_RES, "models")

dir.create(DIR_GEO,  recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PROC, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_META, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TAB,  recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_MOD,  recursive = TRUE, showWarnings = FALSE)
