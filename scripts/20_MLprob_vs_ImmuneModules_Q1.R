# scripts/20_MLprob_vs_ImmuneModules_Q1_UPDATED.R
# Bridge ML outputs -> immune systems biology (CRAN-only)
# What it does:
#  1) Computes immune module scores per sample per dataset (ssGSEA-lite: within-dataset z-score mean)
#  2) Merges ML predicted probabilities (per model) with immune scores (per module)
#  3) Computes Spearman correlation between prob and module score, with BH-FDR across modules (per dataset-model)
#  4) Outputs:
#     - Heatmap of rho (ρ) with stars for each dataset (incl. META pooled)
#     - Scatter panels for BEST model (by pooled AUC) vs key immune modules
#
# Inputs:
#   results/models/STRICT_LODO_multimodel_predictions_long.tsv
#   <DIR_PROC>/<GSE>_expr_gene.rds
#   <DIR_META>/<GSE>_sample_sheet.csv
#
# Outputs:
#   figures/SYSTEMS_MLprob_ImmuneCorr_<tag>/
#     immune_scores_per_sample.tsv
#     pooled_AUC_by_model.tsv
#     MLprob_ImmuneCorr_stats.tsv
#     MLprob_ImmuneCorr_heatmap_<tag>.pdf/png
#     MLprob_vs_ImmuneModules_scatter_bestModel_<BESTMODEL>_<tag>.pdf/png

source("scripts/01_params.R")

# ---------------- Packages (CRAN only) ----------------
options(timeout = 600)
cran_pkgs <- c("data.table","ggplot2","patchwork","stringr","pROC","grid")
for (p in cran_pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(pROC)
  library(grid)
})

options(scipen=999)

# ---------------- Helpers ----------------
DIR_MOD <- if (exists("DIR_MOD")) DIR_MOD else file.path("results","models")
DIR_FIG <- if (exists("DIR_FIG")) DIR_FIG else "figures"
DIR_PROC <- if (exists("DIR_PROC")) DIR_PROC else file.path("results","processed")
DIR_META <- if (exists("DIR_META")) DIR_META else "metadata"

dir.create(DIR_FIG, showWarnings=FALSE, recursive=TRUE)

must_exist <- function(f) if (!file.exists(f)) stop("Missing file: ", f)

safe_tag <- function(x) {
  x <- gsub("<=|≥|≤", "le", x)
  x <- gsub("<", "lt", x)
  x <- gsub(">", "gt", x)
  x <- gsub("=", "eq", x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x
}

p_to_star <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# ---------------- Settings (match your pipeline tag) ----------------
TOPK <- 50
MIN_MODELS <- 3
META_FDR_CUT <- 0.05
tag <- safe_tag(sprintf("TOPK%d_minModels%d_metaFDRle%.2f", TOPK, MIN_MODELS, META_FDR_CUT))

pred_f <- file.path(DIR_MOD, "STRICT_LODO_multimodel_predictions_long.tsv")
must_exist(pred_f)

# ---------------- Load + de-duplicate predictions ----------------
pred0 <- fread(pred_f)
stopifnot(all(c("test_dataset","model","sample_id","y","prob") %in% names(pred0)))

pred0 <- pred0[!is.na(test_dataset) & !is.na(model) & !is.na(sample_id)]
pred0[, y := as.integer(y)]

# Ensure 1 row per dataset-model-sample
pred <- pred0[, .(
  y = max(y, na.rm=TRUE),
  prob = mean(prob, na.rm=TRUE)
), by=.(test_dataset, model, sample_id)]

pred[, class := ifelse(y==1, "HGSC", "Normal")]

DSETS  <- sort(unique(pred$test_dataset))
MODELS <- sort(unique(pred$model))

OUT_DIR <- file.path(DIR_FIG, paste0("SYSTEMS_MLprob_ImmuneCorr_", tag))
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

# Colors
COL_NORMAL <- "#9E9E9E"
COL_HGSC   <- "#E64B35"
COL_POS    <- "#E64B35"
COL_NEG    <- "#00BFC4"

# ---------------- Immune gene sets ----------------
IMMUNE_SETS <- list(
  "IFNG_response"   = c("IFNG","STAT1","IRF1","CXCL9","CXCL10","IDO1","GBP1","GBP2","HLA-DRA","HLA-DRB1"),
  "Cytotoxic_NK"    = c("NKG7","GNLY","PRF1","GZMB","GZMH","CTSW","KLRD1","KLRB1","FCGR3A","TRAC"),
  "Antigen_present" = c("B2M","TAP1","TAP2","HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","CIITA","PSMB9"),
  "T_cell"          = c("CD3D","CD3E","CD3G","TRAC","TRBC1","LCK","ZAP70","LAT","IL7R","CCR7"),
  "B_cell"          = c("MS4A1","CD79A","CD79B","CD74","HLA-DRA","HLA-DRB1","CD37","BANK1","BLK","CD22"),
  "Macrophage"      = c("LST1","TYROBP","AIF1","FCER1G","C1QA","C1QB","C1QC","CSF1R","CTSS","SPI1"),
  "Neutrophil"      = c("S100A8","S100A9","LCN2","MPO","ELANE","FCGR3B","CXCR2","CEACAM8","OLFM4","DEFA4"),
  "Treg"            = c("FOXP3","IL2RA","CTLA4","IKZF2","TIGIT","CCR8","TNFRSF18","ENTPD1","BATF","IRF4"),
  "Checkpoint"      = c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","TIGIT","HAVCR2","TNFRSF9","ICOS","CD80"),
  "Inflammation"    = c("IL1B","IL6","TNF","CXCL8","CCL2","PTGS2","NFKBIA","IRAK3","SOCS3","ICAM1")
)
MODULES <- names(IMMUNE_SETS)

pretty_lab <- c(
  IFNG_response   = "IFN-γ response",
  Cytotoxic_NK    = "Cytotoxic/NK",
  Antigen_present = "Antigen presentation",
  T_cell          = "T cell",
  B_cell          = "B cell",
  Macrophage      = "Macrophage",
  Neutrophil      = "Neutrophil",
  Treg            = "Treg",
  Checkpoint      = "Checkpoint",
  Inflammation    = "Inflammation"
)

# ---------------- Load expression + sheet ----------------
load_expr_and_sheet <- function(ds) {
  expr_path  <- file.path(DIR_PROC, paste0(ds, "_expr_gene.rds"))
  sheet_path <- file.path(DIR_META, paste0(ds, "_sample_sheet.csv"))
  must_exist(expr_path); must_exist(sheet_path)
  
  expr  <- readRDS(expr_path)
  sheet <- fread(sheet_path)
  
  sheet[final_group=="", final_group := NA_character_]
  sheet <- sheet[!is.na(final_group) & final_group %in% c("HGSC","Normal")]
  
  sids <- intersect(colnames(expr), sheet$sample_id)
  if (length(sids) < 6) stop("Too few aligned samples for dataset: ", ds)
  
  sheet <- sheet[match(sids, sample_id)]
  expr  <- expr[, sids, drop=FALSE]
  list(expr=expr, sheet=sheet)
}

zscore_rows <- function(mat) {
  t(apply(mat, 1, function(v){
    m <- mean(v, na.rm=TRUE)
    s <- sd(v, na.rm=TRUE)
    if (!is.finite(s) || s==0) return(rep(0, length(v)))
    (v-m)/s
  }))
}

calc_module_scores <- function(expr, sheet, ds, keep_samples=NULL) {
  if (!is.null(keep_samples)) {
    keep_samples <- intersect(colnames(expr), keep_samples)
    expr <- expr[, keep_samples, drop=FALSE]
    sheet <- sheet[match(keep_samples, sample_id)]
  }
  
  out <- list()
  for (mod in MODULES) {
    gset <- IMMUNE_SETS[[mod]]
    g_avail <- intersect(gset, rownames(expr))
    if (length(g_avail) < 3) next
    
    mat <- expr[g_avail, , drop=FALSE]
    matz <- zscore_rows(mat)
    score <- colMeans(matz, na.rm=TRUE)
    
    out[[mod]] <- data.table(
      test_dataset = ds,
      sample_id = colnames(expr),
      module = mod,
      module_lab = pretty_lab[mod],
      score = as.numeric(score),
      n_genes_used = length(g_avail),
      group_sheet = sheet$final_group
    )
  }
  if (length(out)==0) return(NULL)
  rbindlist(out, use.names=TRUE, fill=TRUE)
}

# ---------------- Compute immune scores ----------------
immune_list <- lapply(DSETS, function(ds) {
  obj <- load_expr_and_sheet(ds)
  keep_s <- unique(pred[test_dataset==ds]$sample_id)
  calc_module_scores(obj$expr, obj$sheet, ds, keep_samples=keep_s)
})
immune0 <- rbindlist(immune_list, use.names=TRUE, fill=TRUE)
if (nrow(immune0)==0) stop("No immune module scores computed. Check gene symbols availability.")

immune <- immune0[!is.na(test_dataset) & !is.na(sample_id) & !is.na(module)]
immune <- immune[, .(
  score = mean(score, na.rm=TRUE),
  n_genes_used = max(n_genes_used, na.rm=TRUE),
  module_lab = first(module_lab),
  group_sheet = first(group_sheet)
), by=.(test_dataset, sample_id, module)]

fwrite(immune, file.path(OUT_DIR, "immune_scores_per_sample.tsv"), sep="\t")

# ---------------- Merge ----------------
dat <- merge(
  pred, immune,
  by=c("test_dataset","sample_id"),
  all.x=TRUE,
  allow.cartesian=TRUE
)
dat <- dat[!is.na(score)]

dat[, dataset := test_dataset]
dat[, module := factor(module, levels=MODULES)]
dat[, module_lab := factor(module_lab, levels=pretty_lab[MODULES])]

dat_meta <- copy(dat)
dat_meta[, dataset := "META (pooled)"]
dat_all <- rbind(dat, dat_meta, use.names=TRUE, fill=TRUE)

dat_all[, dataset := factor(dataset, levels=c(DSETS, "META (pooled)"))]
dat_all[, model := factor(model, levels=MODELS)]
dat_all[, class := factor(class, levels=c("Normal","HGSC"))]

# ---------------- Best model ----------------
auc_tab <- pred[, {
  r <- tryCatch(pROC::roc(y, prob, quiet=TRUE, direction="<"), error=function(e) NULL)
  aucv <- if (is.null(r)) NA_real_ else as.numeric(pROC::auc(r))
  .(AUC=aucv)
}, by=model]
setorder(auc_tab, -AUC)

BEST_MODEL <- as.character(auc_tab[1]$model)
message("Best pooled model by AUC: ", BEST_MODEL)

fwrite(auc_tab, file.path(OUT_DIR, "pooled_AUC_by_model.tsv"), sep="\t")

# ---------------- Correlation stats ----------------
corr_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 8) return(list(rho=NA_real_, p=NA_real_))
  ct <- tryCatch(cor.test(x, y, method="spearman", exact=FALSE), error=function(e) NULL)
  if (is.null(ct)) return(list(rho=NA_real_, p=NA_real_))
  list(rho=as.numeric(ct$estimate), p=ct$p.value)
}

corr <- dat_all[, {
  z <- corr_one(prob, score)
  .(rho=z$rho, p_raw=z$p, n=.N)
}, by=.(dataset, model, module, module_lab)]

corr[, p_adj := p.adjust(p_raw, method="BH"), by=.(dataset, model)]
corr[, star := vapply(p_adj, p_to_star, character(1))]

fwrite(corr, file.path(OUT_DIR, "MLprob_ImmuneCorr_stats.tsv"), sep="\t")

# ---------------- Plot 1: Heatmap ----------------
hm <- copy(corr)
hm[, txt := ifelse(is.na(rho), "", paste0(sprintf("%.2f", rho), "\n", star))]

Lh <- max(abs(hm$rho), na.rm=TRUE)
if (!is.finite(Lh) || Lh==0) Lh <- 1

hm[, module_lab_rev := factor(module_lab, levels=rev(pretty_lab[MODULES]))]

p_hm <- ggplot(hm, aes(x=model, y=module_lab_rev, fill=rho)) +
  geom_tile(color="white", linewidth=0.35) +
  geom_text(
    aes(label = txt),
    color = "black",
    size = 4.5,
    lineheight = 0.95,
    hjust = 0.5,
    vjust = 0.5
  ) +
  scale_fill_gradient2(
    low=COL_NEG, mid="white", high=COL_POS,
    midpoint=0, limits=c(-Lh, Lh), name="Spearman ρ"
  ) +
  facet_wrap(~dataset, ncol=2) +
  theme_bw(base_size=12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(
      size=14, face="bold",
      margin=ggplot2::margin(b=8, unit="pt")
    ),
    strip.text = element_text(size=11, face="bold"),
    axis.title = element_blank(),
    axis.text.x = element_text(size=10, angle=30, hjust=1),
    axis.text.y = element_text(size=10),
    plot.margin = ggplot2::margin(10, 10, 10, 10, unit="pt")
  ) +
  labs(
    title = "Correlation between ML probability and immune module scores",
    subtitle = "Spearman correlation (ρ) with BH-FDR significance per dataset-model (*<0.05, **<0.01, ***<0.001)"
  )

ggsave(
  file.path(OUT_DIR, paste0("MLprob_ImmuneCorr_heatmap_", tag, ".pdf")),
  p_hm, width=12.0, height=9.5, units="in",
  device=cairo_pdf, limitsize=FALSE
)
ggsave(
  file.path(OUT_DIR, paste0("MLprob_ImmuneCorr_heatmap_", tag, ".png")),
  p_hm, width=12.0, height=9.5, units="in",
  dpi=600, bg="white", limitsize=FALSE
)

# ---------------- Plot 2: Scatter panels for best model ----------------
key_modules <- c("Checkpoint","Treg","T_cell","Cytotoxic_NK")
key_labs <- pretty_lab[key_modules]

sc <- dat_all[model == BEST_MODEL & module %in% key_modules]
sc[, module_lab := factor(module_lab, levels = key_labs)]
sc[, dataset := factor(dataset, levels = levels(dat_all$dataset))]

p_sc <- ggplot(sc, aes(x = score, y = prob)) +
  # --- Scatter points (exact size = legend icons)
  geom_point(
    aes(color = class),
    size = 5.6,         # same as legend
    alpha = 0.8,
    stroke = 0.25
  ) +
  # --- Regression line + shaded CI band
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "grey75",
    alpha = 0.45,
    linewidth = 0.9
  ) +
  scale_color_manual(
    values = c("Normal" = COL_NORMAL, "HGSC" = COL_HGSC),
    name = NULL
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 5.6, alpha = 1)
    )
  ) +
  facet_grid(module_lab ~ dataset, scales = "free_x") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      size = 14,
      face = "bold",
      color = "grey30"
    ),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14),
    legend.position = "top",
    legend.key.size = grid::unit(22, "pt"),
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "pt")  # ✅ remove outer margin
  ) +
  labs(
    title = paste0("Best model (pooled AUC): ", BEST_MODEL,
                   " | ML probability vs immune module scores"),
    x = "Immune module score (within-dataset Z-score mean)",
    y = "Predicted probability (ML output)"
  )

# ---- Save ----
ggsave(
  file.path(OUT_DIR, paste0("MLprob_vs_ImmuneModules_scatter_bestModel_", BEST_MODEL, "_", tag, ".pdf")),
  p_sc, width = 16.0, height = 9.5, units = "in",
  device = cairo_pdf, limitsize = FALSE
)
ggsave(
  file.path(OUT_DIR, paste0("MLprob_vs_ImmuneModules_scatter_bestModel_", BEST_MODEL, "_", tag, ".png")),
  p_sc, width = 16.0, height = 9.5, units = "in",
  dpi = 600, bg = "white", limitsize = FALSE
)
