# scripts/09_figures_STRICT_ROC_and_expression_panels.R
# Publication-grade ROC + expression panels (match your provided PDFs)
# Inputs:
#   - results/models/LODO_STRICT_predictions.tsv   (from script 06b)
#   - (optional) results/tables/META_consensus_STRICT.tsv  OR
#               results/models/STRICT_LODO_gene_ranking_aggregated.tsv (from script 11)
# Outputs:
#   - figures/roc_strict_v2/ROC_STRICT_LODO_overlay_AUCfill.pdf
#   - figures/expr_panels_strict_v2/EXPR_UP_top12_panel.pdf
#   - figures/expr_panels_strict_v2/EXPR_DOWN_top12_panel.pdf

source("scripts/01_params.R")

pkgs <- c("data.table","ggplot2","pROC","ggbeeswarm","grid")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pROC)
  library(ggbeeswarm)
  library(grid)
})

# -----------------------------
# paths / dirs
# -----------------------------
DIR_FIG <- if (exists("DIR_FIG")) DIR_FIG else "figures"
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

DIR_ROC  <- file.path(DIR_FIG, "roc_strict_v2")
DIR_EXPR <- file.path(DIR_FIG, "expr_panels_strict_v2")
dir.create(DIR_ROC,  showWarnings=FALSE, recursive=TRUE)
dir.create(DIR_EXPR, showWarnings=FALSE, recursive=TRUE)

DIR_MOD <- if (exists("DIR_MOD")) DIR_MOD else file.path("results","models")
dir.create(DIR_MOD, showWarnings=FALSE, recursive=TRUE)

must_exist <- function(f) if (!file.exists(f)) stop("Missing file: ", f)

# -----------------------------
# helpers
# -----------------------------
p_to_star <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 1e-4) "****"
  else if (p < 1e-3) "***"
  else if (p < 1e-2) "**"
  else if (p < 0.05) "*"
  else "ns"
}

zscore <- function(x) {
  m <- mean(x, na.rm=TRUE)
  s <- sd(x, na.rm=TRUE)
  if (!is.finite(s) || s == 0) return(x*0)
  (x - m) / s
}

load_expr_and_sheet <- function(ds) {
  expr_path  <- file.path(DIR_PROC, paste0(ds, "_expr_gene.rds"))
  sheet_path <- file.path(DIR_META, paste0(ds, "_sample_sheet.csv"))
  must_exist(expr_path); must_exist(sheet_path)
  
  expr  <- readRDS(expr_path)           # genes x samples
  sheet <- fread(sheet_path)
  
  if (!("sample_id" %in% names(sheet))) stop("sample_id missing in ", sheet_path)
  if (!("final_group" %in% names(sheet))) stop("final_group missing in ", sheet_path)
  
  sheet$final_group[sheet$final_group==""] <- NA
  sheet <- sheet[!is.na(final_group) & final_group %in% c("HGSC","Normal")]
  
  sids <- intersect(colnames(expr), sheet$sample_id)
  if (length(sids) < 6) stop("Too few matched samples in ", ds)
  
  sheet <- sheet[match(sids, sheet$sample_id)]
  expr  <- expr[, sids, drop=FALSE]
  
  list(expr=expr, sheet=sheet)
}

# -----------------------------
# 1) ROC overlay (outline + linetype + micro-offset) + AUC labels + light fill
# -----------------------------
pred_f <- file.path(DIR_MOD, "LODO_STRICT_predictions.tsv")
must_exist(pred_f)
pred <- fread(pred_f)
pred$y <- as.integer(pred$y)
pred <- pred[!is.na(y) & !is.na(prob)]

dsets <- sort(unique(pred$test_dataset))

# compute AUC per dataset
auc_tab <- rbindlist(lapply(dsets, function(ds) {
  dd <- pred[test_dataset==ds]
  r  <- pROC::roc(dd$y, dd$prob, quiet=TRUE, direction="<")
  data.table(dataset=ds, AUC=as.numeric(pROC::auc(r)))
}))

label_map <- setNames(
  sprintf("%s (AUC=%.3f)", auc_tab$dataset, auc_tab$AUC),
  auc_tab$dataset
)

roc_df_from_pred <- function(dd) {
  r <- pROC::roc(dd$y, dd$prob, quiet=TRUE, direction="<")
  coords_all <- pROC::coords(r, "all", ret=c("specificity","sensitivity"), transpose=FALSE)
  data.table(
    FPR = 1 - coords_all$specificity,
    TPR = coords_all$sensitivity
  )
}

roc_list <- list()
for (i in seq_along(dsets)) {
  ds <- dsets[i]
  dd <- pred[test_dataset==ds]
  rdf <- roc_df_from_pred(dd)
  
  # micro-offset only for visibility when curves overlap
  eps <- 0.001
  rdf[, FPR := pmin(pmax(FPR + (i-1)*eps, 0), 1)]
  
  rdf[, dataset := ds]
  rdf[, dataset_lbl := label_map[ds]]
  roc_list[[ds]] <- rdf
}
rocdf <- rbindlist(roc_list)
rocdf[, dataset_lbl := factor(dataset_lbl, levels=label_map[dsets])]

p_roc <- ggplot(rocdf, aes(x=FPR, y=TPR, group=dataset_lbl, color=dataset_lbl, linetype=dataset_lbl)) +
  # light fill under curve
  geom_ribbon(aes(ymin=0, ymax=TPR, fill=dataset_lbl), alpha=0.08, color=NA, show.legend=FALSE) +
  # diagonal
  geom_abline(intercept=0, slope=1, linewidth=0.7, linetype="dashed", color="grey55") +
  # outline trick (black behind)
  geom_path(linewidth=3.0, color="black", alpha=0.55, lineend="round") +
  geom_path(linewidth=1.6, lineend="round") +
  coord_equal() +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color="grey90", linewidth=0.35),
    legend.position  = c(0.78, 0.20),
    legend.title     = element_text(size=11),
    legend.text      = element_text(size=10),
    legend.key.width = unit(1.2, "cm"),
    plot.title       = element_text(size=18, face="bold", hjust=0.5),
    axis.title       = element_text(size=14)
  ) +
  labs(
    title="STRICT LODO ROC (Overlay: outline + linetype + micro-offset)",
    x="False Positive Rate",
    y="True Positive Rate",
    color="Test dataset",
    linetype="Test dataset"
  )

roc_pdf <- file.path(DIR_ROC, "ROC_STRICT_LODO_overlay_AUCfill.pdf")
ggsave(roc_pdf, p_roc, width=7.6, height=7.0)
message("Saved: ", roc_pdf)

# -----------------------------
# 2) Choose TOP genes (prefer META consensus; fallback to aggregated ranking from script 11)
# -----------------------------
meta_candidates <- c(
  file.path(DIR_TAB, "META_consensus_STRICT.tsv"),
  file.path(DIR_TAB, "META_consensus_RELAXED.tsv"),
  file.path(DIR_MOD, "STRICT_LODO_gene_ranking_aggregated.tsv")
)
meta_f <- meta_candidates[file.exists(meta_candidates)][1]
if (is.na(meta_f)) stop("No meta/ranking table found. Expected one of: ", paste(meta_candidates, collapse=" | "))
meta <- fread(meta_f)

# standardize columns
if (!("gene" %in% names(meta))) stop("Meta/ranking table lacks 'gene' column: ", meta_f)

# find logFC column
logfc_col <- c("meta_logFC","logFC","Meta_logFC","META_logFC")[c("meta_logFC","logFC","Meta_logFC","META_logFC") %in% names(meta)][1]
# find FDR column
fdr_col <- c("meta_fdr","meta_FDR","FDR","padj","adj.P.Val","fdr")[c("meta_fdr","meta_FDR","FDR","padj","adj.P.Val","fdr") %in% names(meta)][1]

if (is.na(logfc_col) || is.na(fdr_col)) {
  stop("Could not detect meta_logFC/meta_fdr columns in: ", meta_f,
       "\nFound columns: ", paste(names(meta), collapse=", "))
}

meta[, meta_logFC := as.numeric(get(logfc_col))]
meta[, meta_fdr   := as.numeric(get(fdr_col))]
meta <- meta[is.finite(meta_logFC) & is.finite(meta_fdr)]

# direction
meta[, direction := ifelse(meta_logFC > 0, "UP", "DOWN")]
meta[, abs_logFC := abs(meta_logFC)]

# select top genes
TOP_UP <- 12
TOP_DN <- 12
FDR_CUT <- 0.05

meta_sig <- meta[meta_fdr <= FDR_CUT]

up <- meta_sig[direction=="UP"]
dn <- meta_sig[direction=="DOWN"]

# robust ordering (NO get()/abs() inside setorder)
setorder(up, meta_fdr, -abs_logFC)
setorder(dn, meta_fdr, -abs_logFC)

up_genes <- head(up$gene, TOP_UP)
dn_genes <- head(dn$gene, TOP_DN)

message("[INFO] Gene source: ", basename(meta_f))
message("[INFO] Selected UP genes: ", length(up_genes), " | DOWN genes: ", length(dn_genes))

# datasets from prediction file (your 3 GEO sets)
datasets <- sort(unique(pred$test_dataset))

# -----------------------------
# 3) Build expression panel: gene rows x dataset columns (3 GEO + META)
# -----------------------------
make_expr_long <- function(genes, datasets) {
  rows <- list()
  
  for (ds in datasets) {
    obj <- load_expr_and_sheet(ds)
    expr <- obj$expr
    sheet <- obj$sheet
    
    keep_genes <- intersect(genes, rownames(expr))
    if (length(keep_genes)==0) next
    
    mat <- expr[keep_genes, , drop=FALSE]
    # z-score within dataset (per gene)
    zmat <- t(apply(mat, 1, zscore))
    
    dd <- data.table(
      gene  = rep(rownames(zmat), times=ncol(zmat)),
      dataset = ds,
      sample  = rep(colnames(zmat), each=nrow(zmat)),
      group = rep(sheet$final_group, each=nrow(zmat)),
      value = as.numeric(as.vector(zmat))
    )
    rows[[ds]] <- dd
  }
  
  df <- rbindlist(rows, fill=TRUE)
  if (nrow(df)==0) return(NULL)
  
  # META pool (already within-dataset z-scored)
  meta_df <- copy(df)
  meta_df[, dataset := "META"]
  rbind(df, meta_df)
}

annot_stars <- function(df) {
  # wilcox per (gene, dataset): Normal vs HGSC
  pv <- df[, .(
    p = if (length(unique(group)) < 2) NA_real_ else wilcox.test(value ~ group)$p.value
  ), by=.(gene, dataset)]
  pv[, star := vapply(p, p_to_star, character(1))]
  
  y_max <- df[, .(ymax=max(value, na.rm=TRUE)), by=.(gene, dataset)]
  ann <- merge(pv, y_max, by=c("gene","dataset"), all.x=TRUE)
  ann[, y := ymax + 0.35]
  ann
}

plot_expr_panel <- function(genes, datasets, title, out_pdf) {
  df <- make_expr_long(genes, datasets)
  if (is.null(df)) {
    message("[SKIP] No expression rows for this gene set.")
    return(NULL)
  }
  
  # ordering
  ds_levels <- c(datasets, "META")
  df[, dataset := factor(dataset, levels=ds_levels)]
  df[, group := factor(group, levels=c("Normal","HGSC"))]
  
  # order genes by meta significance (if available)
  ord_tbl <- meta[gene %in% genes, .(gene, meta_fdr, abs_logFC)]
  setorder(ord_tbl, meta_fdr, -abs_logFC)
  gene_levels <- ord_tbl$gene
  df[, gene := factor(gene, levels=rev(gene_levels))]  # top at top like your PDF
  
  ann <- annot_stars(df)
  ann[, dataset := factor(dataset, levels=ds_levels)]
  ann[, gene := factor(gene, levels=levels(df$gene))]
  
  # short x labels to avoid overlap
  x_labs <- c("Normal"="N", "HGSC"="HGSC")
  
  # color palette (as in your PDFs)
  pal_fill <- c("Normal"="#BDBDBD", "HGSC"="#4C72B0")
  
  ng <- length(unique(df$gene))
  h  <- max(5.2, 0.55*ng + 2.0)
  
  p <- ggplot(df, aes(x=group, y=value, fill=group)) +
    geom_violin(trim=TRUE, color=NA, alpha=0.55, scale="width") +
    geom_boxplot(width=0.16, outlier.shape=NA, color=NA, alpha=0.75) +
    ggbeeswarm::geom_quasirandom(size=0.65, color="black", width=0.18, alpha=0.80) +
    facet_grid(gene ~ dataset) +
    geom_text(data=ann, aes(x=1.5, y=y, label=star), inherit.aes=FALSE, size=2.8) +
    scale_x_discrete(labels=x_labs) +
    scale_fill_manual(values=pal_fill) +
    theme_bw(base_size=10) +
    theme(
      legend.position="none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color="grey92", linewidth=0.30),
      strip.background = element_rect(fill="white", color="grey85", linewidth=0.25),
      strip.text.x = element_text(size=10),
      strip.text.y = element_text(size=9, face="italic"),
      axis.text.x  = element_text(size=8),
      axis.text.y  = element_text(size=8),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=10),
      panel.spacing = unit(0.18, "lines"),
      plot.title = element_text(size=14, face="bold", hjust=0.5),
      plot.margin = margin(8, 10, 8, 10),
      panel.border = element_rect(color="grey70", linewidth=0.25)
    ) +
    labs(title=title, y="zExpr")
  
  ggsave(out_pdf, p, width=13.0, height=h)
  message("Saved: ", out_pdf)
  p
}

# -----------------------------
# 4) Save UP/DOWN panels (match your PDFs)
# -----------------------------
if (length(up_genes) > 0) {
  plot_expr_panel(
    genes=up_genes,
    datasets=datasets,
    title=sprintf("Top UP genes (metaFDR ≤ %.2f) | z-score within dataset", FDR_CUT),
    out_pdf=file.path(DIR_EXPR, "EXPR_UP_top12_panel.pdf")
  )
} else {
  message("[WARN] No UP genes passed metaFDR cutoff.")
}

if (length(dn_genes) > 0) {
  plot_expr_panel(
    genes=dn_genes,
    datasets=datasets,
    title=sprintf("Top DOWN genes (metaFDR ≤ %.2f) | z-score within dataset", FDR_CUT),
    out_pdf=file.path(DIR_EXPR, "EXPR_DOWN_top12_panel.pdf")
  )
} else {
  message("[WARN] No DOWN genes passed metaFDR cutoff.")
}

message("\nDONE (script 09).")
