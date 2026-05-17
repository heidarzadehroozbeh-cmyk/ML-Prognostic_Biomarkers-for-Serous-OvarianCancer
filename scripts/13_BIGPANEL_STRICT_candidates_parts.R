# scripts/13_BIGPANEL_STRICT_candidates_parts.R
# Export each component of BIGPANEL as separate clean figures (Q1-ready, no overlap).
# Outputs: PDF (vector, cairo) + PNG (600 dpi) per panel.

source("scripts/01_params.R")

pkgs <- c("data.table","ggplot2","pROC","stringr","Matrix","glmnet","matrixStats","grid")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
if (!requireNamespace("ggVennDiagram", quietly=TRUE)) install.packages("ggVennDiagram")

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(pROC); library(stringr)
  library(Matrix); library(glmnet); library(matrixStats); library(ggVennDiagram); library(grid)
})

options(scipen=999)

DIR_MOD <- if (exists("DIR_MOD")) DIR_MOD else file.path("results","models")
DIR_FIG <- if (exists("DIR_FIG")) DIR_FIG else "figures"
dir.create(DIR_MOD, showWarnings=FALSE, recursive=TRUE)
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

TOPK <- 50
MIN_MODELS <- 3
META_FDR_CUT <- 0.05
tag <- safe_tag(sprintf("TOPK%d_minModels%d_metaFDRle%.2f", TOPK, MIN_MODELS, META_FDR_CUT))

OUT_DIR <- file.path(DIR_FIG, paste0("BIGPANEL_PARTS_", tag))
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

pred_long_f <- file.path(DIR_MOD, "STRICT_LODO_multimodel_predictions_long.tsv")
imp_long_f  <- file.path(DIR_MOD, "STRICT_LODO_multimodel_importance_long.tsv")
rank_ag_f   <- file.path(DIR_MOD, "STRICT_LODO_gene_ranking_aggregated.tsv")

must_exist(pred_long_f); must_exist(imp_long_f); must_exist(rank_ag_f)

pred_long <- fread(pred_long_f)
imp_long  <- fread(imp_long_f)
rank_ag   <- fread(rank_ag_f)

cand_fin_f <- file.path(DIR_MOD, paste0("STRICT_candidates_FINAL_", tag, ".tsv"))
topk_f     <- file.path(DIR_MOD, paste0("STRICT_candidates_topk_by_model_", tag, ".tsv"))

final_list <- if (file.exists(cand_fin_f)) fread(cand_fin_f) else NULL
topk_by_model <- if (file.exists(topk_f)) fread(topk_f) else NULL

if (is.null(topk_by_model)) {
  imp_ag <- imp_long[, .(mean_importance = mean(as.numeric(importance), na.rm=TRUE)), by=.(model, gene)]
  imp_ag[, rank_in_model := frank(-mean_importance, ties.method="average"), by=model]
  topk_by_model <- imp_ag[rank_in_model <= TOPK]
}

models <- sort(unique(pred_long$model))
dsets  <- sort(unique(pred_long$test_dataset))

base_ds <- c("#E64B35","#00A087","#3C5488","#7E6148","#4DBBD5")
cols_ds <- setNames(base_ds[seq_along(dsets)], dsets)

base_m <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")
model_cols <- setNames(base_m[seq_along(models)], models)

save_dual <- function(p, stem, w=7, h=5, dpi=600) {
  pdf_f <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  png_f <- file.path(OUT_DIR, paste0(stem, ".png"))
  ggsave(pdf_f, p, width=w, height=h, units="in", device=cairo_pdf, limitsize=FALSE)
  ggsave(png_f, p, width=w, height=h, units="in", dpi=dpi, bg="white", limitsize=FALSE)
}

roc_df_one <- function(dd, offset_idx=0, eps=0.001) {
  r <- pROC::roc(dd$y, dd$prob, quiet=TRUE, direction="<")
  coords_all <- pROC::coords(r, "all", ret=c("specificity","sensitivity"), transpose=FALSE)
  out <- data.table(
    FPR = 1 - coords_all$specificity,
    TPR = coords_all$sensitivity
  )
  out[, FPR := pmin(pmax(FPR + offset_idx*eps, 0), 1)]
  out
}

plot_roc_model <- function(model_name) {
  
  ddm <- pred_long[model==model_name]
  if (nrow(ddm)==0) return(NULL)
  
  ddm[, y := as.integer(y)]
  
  roc_list <- list()
  auc_tab  <- list()
  
  for (i in seq_along(dsets)) {
    
    ds <- dsets[i]
    dd <- ddm[test_dataset==ds]
    
    if (nrow(dd)==0) next
    
    rdf <- roc_df_one(dd, offset_idx=i-1, eps=0.001)
    rdf[, dataset := ds]
    roc_list[[ds]] <- rdf
    
    r <- pROC::roc(dd$y, dd$prob, quiet=TRUE, direction="<")
    auc_tab[[ds]] <- data.table(model=model_name, dataset=ds, AUC=as.numeric(pROC::auc(r)))
  }
  
  rocdf <- rbindlist(roc_list, use.names=TRUE, fill=TRUE)
  aucdf <- rbindlist(auc_tab,  use.names=TRUE, fill=TRUE)
  
  if (nrow(rocdf)==0) return(NULL)
  
  aucdf[, dataset_lab := sprintf("%s (AUC=%.3f)", dataset, AUC)]
  rocdf <- merge(rocdf, aucdf[, .(dataset, dataset_lab)], by="dataset", all.x=TRUE)
  
  p <- ggplot(rocdf, aes(FPR, TPR, group=dataset, color=dataset, linetype=dataset)) +
    geom_abline(intercept=0, slope=1, linewidth=0.35, linetype="dashed", color="grey55") +
    geom_ribbon(aes(ymin=0, ymax=TPR, fill=dataset), alpha=0.08, color=NA) +
    geom_path(linewidth=1.15) +
    coord_equal(clip="off") +
    scale_color_manual(values=cols_ds, breaks=dsets,
                       labels=aucdf$dataset_lab[match(dsets, aucdf$dataset)]) +
    scale_fill_manual(values=cols_ds, guide="none") +
    theme_bw(base_size=11) +
    theme(
      plot.title = element_text(size=13, face="bold",
                                margin=ggplot2::margin(b=6, unit="pt")),
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.text  = element_text(size=9),
      legend.key.width = unit(1.8, "lines"),
      panel.grid.minor = element_blank(),
      plot.margin = ggplot2::margin(10,10,10,10,"pt")
    ) +
    labs(title=paste0("ROC – ", model_name),
         x="False Positive Rate", y="True Positive Rate", color=NULL, linetype=NULL)
  
  list(plot=p, auc=aucdf)
}

auc_all <- list()
for (m in models) {
  tmp <- plot_roc_model(m)
  if (is.null(tmp)) next
  save_dual(tmp$plot, stem=paste0("A_ROC_", make.names(m)), w=7.2, h=6.2)
  auc_all[[m]] <- tmp$auc
}

auc_tbl <- rbindlist(auc_all, use.names=TRUE, fill=TRUE)
fwrite(auc_tbl, file.path(OUT_DIR, "ROC_AUC_table.tsv"), sep="\t")

message("DONE. Saved individual panels to: ", OUT_DIR)
