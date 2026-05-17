# scripts/12_select_candidate_biomarkers_STRICT_LODO.R
# Candidate biomarker selection using CLASSIC ML models (ElasticNet / SVM / RF / XGBoost / GBM)
# Inputs (from script 11):
#  - results/models/STRICT_LODO_multimodel_importance_long.tsv
#  - results/models/STRICT_LODO_gene_ranking_aggregated.tsv   (preferred)
# Outputs:
#  - results/models/STRICT_candidates_topk_by_model_*.tsv
#  - results/models/STRICT_candidates_all_support_*.tsv
#  - results/models/STRICT_candidates_FINAL_*.tsv
#  - figures/candidates_viz/Venn_TOPK_by_model_*.pdf (+png)
#  - figures/candidates_viz/SupportHeatmap_FINAL_*.pdf
#  - figures/candidates_viz/SupportBar_FINAL_*.pdf

source("scripts/01_params.R")

pkgs <- c("data.table","ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------
# Parameters (edit here)
# -----------------------------
TOPK_PER_MODEL <- 50        # top genes per model (based on aggregated importance)
MIN_MODELS     <- 3         # gene must appear in >= this many model-topK sets
META_FDR_CUT   <- 0.05      # meta significance filter
ALLOW_NA_META  <- FALSE     # if FALSE: genes with NA meta_fdr are excluded from FINAL

# Which models (must match "model" field in importance_long)
MODELS <- c("ElasticNet","SVM","RandomForest","XGBoost","GBM")

# -----------------------------
# Paths
# -----------------------------
DIR_MOD <- if (exists("DIR_MOD")) DIR_MOD else file.path("results","models")
DIR_FIG <- if (exists("DIR_FIG")) DIR_FIG else "figures"

dir.create(DIR_MOD, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

DIR_VIZ <- file.path(DIR_FIG, "candidates_viz")
dir.create(DIR_VIZ, showWarnings = FALSE, recursive = TRUE)

f_imp  <- file.path(DIR_MOD, "STRICT_LODO_multimodel_importance_long.tsv")
f_rank <- file.path(DIR_MOD, "STRICT_LODO_gene_ranking_aggregated.tsv")

must_exist <- function(f) if (!file.exists(f)) stop("Missing file: ", f)

must_exist(f_imp)
must_exist(f_rank)

# -----------------------------
# Safe filename (Windows-proof)
#  - replaces illegal characters: <>:"/\|?*
# -----------------------------
safe_fname <- function(x) {
  x <- gsub("[<>:\"/\\\\|?*]", "_", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("\\.$", "", x)
  x
}

tag_raw <- sprintf("TOPK%d_minModels%d_metaFDRle%.2g", TOPK_PER_MODEL, MIN_MODELS, META_FDR_CUT)
tag <- safe_fname(tag_raw)

f_topk <- file.path(DIR_MOD, paste0("STRICT_candidates_topk_by_model_", tag, ".tsv"))
f_all  <- file.path(DIR_MOD, paste0("STRICT_candidates_all_support_", tag, ".tsv"))
f_fin  <- file.path(DIR_MOD, paste0("STRICT_candidates_FINAL_", tag, ".tsv"))

# -----------------------------
# Read inputs
# -----------------------------
imp <- fread(f_imp)
rank <- fread(f_rank)

# Basic checks
need_cols_imp <- c("test_dataset","model","gene","importance")
if (!all(need_cols_imp %in% names(imp))) {
  stop("importance_long missing required columns. Need: ", paste(need_cols_imp, collapse=", "),
       " | Found: ", paste(names(imp), collapse=", "))
}
if (!("gene" %in% names(rank))) stop("ranking file missing 'gene' column.")

# Harmonize model names
imp$model <- as.character(imp$model)
imp$gene  <- as.character(imp$gene)

# Keep only requested models if present
present_models <- intersect(MODELS, unique(imp$model))
if (length(present_models) < 2) stop("Not enough models found in importance_long. Found: ", paste(unique(imp$model), collapse=", "))

MODELS <- present_models

# Meta columns (optional but expected)
# ranking file from script 11 typically has: meta_logFC, meta_fdr, direction
meta_cols <- c("meta_logFC","meta_fdr","direction")
for (cc in meta_cols) if (!(cc %in% names(rank))) rank[, (cc) := NA]

# Ensure numeric
rank[, meta_logFC := as.numeric(meta_logFC)]
rank[, meta_fdr   := as.numeric(meta_fdr)]

# -----------------------------
# 1) Build TOPK sets per model (aggregate importance across folds)
# -----------------------------
topk_by_model_list <- list()

for (m in MODELS) {
  dd <- imp[model == m & is.finite(importance)]
  if (nrow(dd) == 0) next
  
  agg <- dd[, .(
    mean_importance = mean(importance, na.rm=TRUE),
    median_importance = median(importance, na.rm=TRUE),
    n_entries = .N
  ), by=gene]
  
  setorder(agg, -mean_importance, -median_importance, -n_entries, gene)
  
  k <- min(TOPK_PER_MODEL, nrow(agg))
  agg <- agg[1:k]
  agg[, model := m]
  agg[, rank_in_model := .I]
  
  topk_by_model_list[[m]] <- agg
}

topk_by_model <- rbindlist(topk_by_model_list, fill=TRUE)
if (nrow(topk_by_model) == 0) stop("topk_by_model is empty. Check importance_long content.")

# Attach meta info
topk_by_model <- merge(
  topk_by_model,
  unique(rank[, .(gene, meta_logFC, meta_fdr, direction)]),
  by="gene",
  all.x=TRUE
)

# Save
fwrite(topk_by_model, f_topk, sep="\t")

# -----------------------------
# 2) Model-support table (how many models select the gene in TOPK)
# -----------------------------
support <- topk_by_model[, .(
  n_models = uniqueN(model),
  models   = paste(sort(unique(model)), collapse=";"),
  best_rank = min(rank_in_model, na.rm=TRUE),
  mean_importance = mean(mean_importance, na.rm=TRUE),
  meta_logFC = meta_logFC[which.min(is.na(meta_logFC))][1],
  meta_fdr   = meta_fdr[which.min(is.na(meta_fdr))][1],
  direction  = direction[which.min(is.na(direction))][1]
), by=gene]

# Order for readability
# data.table cannot sort on expression abs(meta_logFC) directly → create helper column
support[, abs_meta_logFC := abs(meta_logFC)]
setorder(support, -n_models, meta_fdr, -abs_meta_logFC, best_rank, gene)

fwrite(support, f_all, sep="\t")

# -----------------------------
# 3) Final candidate list (filters)
# -----------------------------
final_list <- copy(support)
final_list <- final_list[n_models >= MIN_MODELS]

if (!ALLOW_NA_META) {
  final_list <- final_list[!is.na(meta_fdr)]
}
final_list <- final_list[is.na(meta_fdr) | meta_fdr <= META_FDR_CUT]

# Re-order final
final_list[, abs_meta_logFC := abs(meta_logFC)]
setorder(final_list, -n_models, meta_fdr, -abs_meta_logFC, best_rank, gene)

# remove helper
final_list[, abs_meta_logFC := NULL]
support[,    abs_meta_logFC := NULL]

fwrite(final_list, f_fin, sep="\t")

message("\nSaved:")
message(" - ", f_topk)
message(" - ", f_all)
message(" - ", f_fin)

message("\nFinal candidates: n=", nrow(final_list))
if (nrow(final_list) > 0) {
  print(head(final_list, 30))
}

# -----------------------------
# 4) Venn diagram (nice) for TOPK sets by model
# -----------------------------
# Build gene sets per model
gene_sets <- lapply(MODELS, function(m) unique(topk_by_model[model == m]$gene))
names(gene_sets) <- MODELS

# Prefer ggVennDiagram; fallback to VennDiagram
venn_pdf <- file.path(DIR_VIZ, paste0("Venn_TOPK_by_model_", tag, ".pdf"))
venn_png <- file.path(DIR_VIZ, paste0("Venn_TOPK_by_model_", tag, ".png"))

venn_done <- FALSE

if (requireNamespace("ggVennDiagram", quietly=TRUE)) {
  suppressPackageStartupMessages(library(ggVennDiagram))
  
  pvenn <- ggVennDiagram::ggVennDiagram(
    gene_sets,
    label = "count",
    label_alpha = 0
  ) +
    ggtitle(sprintf("Top-%d genes per model (counts shown)", TOPK_PER_MODEL)) +
    theme(plot.title = element_text(size=12, face="bold"))
  
  ggsave(venn_pdf, pvenn, width=8.0, height=7.0)
  ggsave(venn_png, pvenn, width=8.0, height=7.0, dpi=300)
  venn_done <- TRUE
} else {
  # try install
  try({
    install.packages("ggVennDiagram")
  }, silent=TRUE)
  
  if (requireNamespace("ggVennDiagram", quietly=TRUE)) {
    suppressPackageStartupMessages(library(ggVennDiagram))
    pvenn <- ggVennDiagram::ggVennDiagram(gene_sets, label = "count", label_alpha = 0) +
      ggtitle(sprintf("Top-%d genes per model (counts shown)", TOPK_PER_MODEL)) +
      theme(plot.title = element_text(size=12, face="bold"))
    ggsave(venn_pdf, pvenn, width=8.0, height=7.0)
    ggsave(venn_png, pvenn, width=8.0, height=7.0, dpi=300)
    venn_done <- TRUE
  }
}

if (!venn_done) {
  if (!requireNamespace("VennDiagram", quietly=TRUE)) install.packages("VennDiagram")
  suppressPackageStartupMessages(library(VennDiagram))
  
  # VennDiagram draws to a device; produce pdf
  pdf(venn_pdf, width=8.5, height=7.5)
  grid::grid.newpage()
  VennDiagram::venn.diagram(
    x = gene_sets,
    filename = NULL,
    main = sprintf("Top-%d genes per model (counts shown)", TOPK_PER_MODEL),
    main.cex = 1.1,
    cat.cex = 0.9,
    cex = 0.9
  ) |> grid::grid.draw()
  dev.off()
  
  # png too
  png(venn_png, width=2400, height=2000, res=300)
  grid::grid.newpage()
  VennDiagram::venn.diagram(
    x = gene_sets,
    filename = NULL,
    main = sprintf("Top-%d genes per model (counts shown)", TOPK_PER_MODEL),
    main.cex = 1.1,
    cat.cex = 0.9,
    cex = 0.9
  ) |> grid::grid.draw()
  dev.off()
}

message("Saved Venn: ", venn_pdf)
message("Saved Venn: ", venn_png)

# -----------------------------
# 5) Support heatmap (FINAL) + barplot
# -----------------------------
if (nrow(final_list) > 0) {
  
  # Build matrix (gene x model) from TOPK membership
  topk_flag <- unique(topk_by_model[, .(model, gene)])
  topk_flag[, selected := 1L]
  
  # restrict to final genes for cleaner figure
  genes_plot <- final_list$gene
  mm <- CJ(model = MODELS, gene = genes_plot)
  mm <- merge(mm, topk_flag, by=c("model","gene"), all.x=TRUE)
  mm[is.na(selected), selected := 0L]
  
  # join meta for ordering (UP then DOWN; higher support first)
  meta_plot <- unique(final_list[, .(gene, n_models, meta_fdr, meta_logFC, direction)])
  meta_plot[, abs_meta_logFC := abs(meta_logFC)]
  # order genes: support desc, meta_fdr asc, abs(meta_logFC) desc
  setorder(meta_plot, -n_models, meta_fdr, -abs_meta_logFC, gene)
  
  mm <- merge(mm, meta_plot, by="gene", all.x=TRUE)
  
  # Heatmap
  mm$gene <- factor(mm$gene, levels=meta_plot$gene)
  mm$model <- factor(mm$model, levels=MODELS)
  
  p_hm <- ggplot(mm, aes(x=model, y=gene, fill=factor(selected))) +
    geom_tile(color="grey85", linewidth=0.25) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size=7),
      legend.position = "none"
    ) +
    labs(
      title = sprintf("FINAL candidates: membership in Top-%d per model", TOPK_PER_MODEL),
      x = NULL, y = NULL
    )
  
  hm_pdf <- file.path(DIR_VIZ, paste0("SupportHeatmap_FINAL_", tag, ".pdf"))
  ggsave(hm_pdf, p_hm, width=6.2, height=max(4.5, 0.18*length(genes_plot)))
  message("Saved: ", hm_pdf)
  
  # Support bar
  p_bar <- ggplot(meta_plot, aes(x=reorder(gene, n_models), y=n_models)) +
    geom_col() +
    coord_flip() +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank()) +
    labs(title="FINAL candidates: model-support count", x=NULL, y="#Models (TopK membership)")
  
  bar_pdf <- file.path(DIR_VIZ, paste0("SupportBar_FINAL_", tag, ".pdf"))
  ggsave(bar_pdf, p_bar, width=6.2, height=max(4.5, 0.18*length(genes_plot)))
  message("Saved: ", bar_pdf)
}

message("\nDone script 12.")
