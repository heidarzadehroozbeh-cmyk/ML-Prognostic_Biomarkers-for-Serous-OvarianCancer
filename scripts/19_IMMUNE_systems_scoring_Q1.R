# scripts/19_IMMUNE_systems_scoring_Q1.R
# Immune systems biology (CRAN-only)

source("scripts/01_params.R")

options(timeout = 600)
cran_pkgs <- c("data.table","ggplot2","patchwork","stringr")
for (p in cran_pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(stringr)
})

options(scipen=999)
set.seed(1)

DIR_FIG <- if (exists("DIR_FIG")) DIR_FIG else "figures"
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

DSETS <- c("GSE14407","GSE38666","GSE52037")

OUT_DIR <- file.path(DIR_FIG, paste0("SYSTEMS_ImmuneModules_", tag))
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

COL_NORMAL <- "#9E9E9E"
COL_UPWARM <- "#E64B35"
COL_DOWNTR <- "#00BFC4"

p_to_star <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

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

load_expr_and_sheet <- function(ds) {
  expr_path  <- file.path(DIR_PROC, paste0(ds, "_expr_gene.rds"))
  sheet_path <- file.path(DIR_META, paste0(ds, "_sample_sheet.csv"))
  must_exist(expr_path); must_exist(sheet_path)
  
  expr  <- readRDS(expr_path)
  sheet <- fread(sheet_path)
  sheet[final_group=="", final_group := NA_character_]
  sheet <- sheet[!is.na(final_group) & final_group %in% c("HGSC","Normal")]
  
  sids <- intersect(colnames(expr), sheet$sample_id)
  sheet <- sheet[match(sids, sample_id)]
  expr  <- expr[, sids, drop=FALSE]
  
  list(expr=expr, sheet=sheet)
}

zscore_rows <- function(mat) {
  t(apply(mat, 1, function(v){
    m <- mean(v, na.rm=TRUE); s <- sd(v, na.rm=TRUE)
    if (!is.finite(s) || s==0) return(rep(0, length(v)))
    (v-m)/s
  }))
}

calc_module_scores <- function(expr, sheet, ds) {
  out <- list()
  for (mod in MODULES) {
    g <- intersect(IMMUNE_SETS[[mod]], rownames(expr))
    if (length(g) < 3) next
    matz <- zscore_rows(expr[g,,drop=FALSE])
    score <- colMeans(matz)
    out[[mod]] <- data.table(
      sample_id = colnames(expr),
      dataset=ds,
      group=sheet$final_group,
      module=mod,
      module_lab=pretty_lab[mod],
      score,
      n_genes_used=length(g)
    )
  }
  rbindlist(out)
}

scores <- rbindlist(lapply(DSETS, function(ds){
  obj <- load_expr_and_sheet(ds)
  calc_module_scores(obj$expr, obj$sheet, ds)
}))

scores_meta <- copy(scores)[, dataset := "META (pooled)"]
scores_all <- rbind(scores, scores_meta)

scores_all[, group := factor(group, levels=c("Normal","HGSC"))]
scores_all[, module_lab := factor(module_lab, levels=pretty_lab[MODULES])]
scores_all[, dataset := factor(dataset, levels=c(DSETS,"META (pooled)"))]

boot_ci_median_diff <- function(xe, xn, B=300) {
  xe <- xe[is.finite(xe)]
  xn <- xn[is.finite(xn)]
  diffs <- replicate(B,{
    median(sample(xe, length(xe), TRUE)) -
      median(sample(xn, length(xn), TRUE))
  })
  quantile(diffs, c(0.025,0.975))
}

stats_one_ds <- function(ds) {
  dd <- scores_all[dataset==ds]
  st <- dd[, .(
    n_norm=sum(group=="Normal"),
    n_HGSC=sum(group=="HGSC"),
    med_norm=median(score[group=="Normal"]),
    med_HGSC=median(score[group=="HGSC"]),
    delta=median(score[group=="HGSC"]) - median(score[group=="Normal"]),
    p_raw=suppressWarnings(wilcox.test(score~group)$p.value)
  ), by=.(module, module_lab)]
  
  st[, p_adj := p.adjust(p_raw,"BH")]
  st[, star := vapply(p_adj,p_to_star,"")]
  
  ci <- lapply(st$module, function(m){
    xe <- dd[module==m & group=="HGSC"]$score
    xn <- dd[module==m & group=="Normal"]$score
    boot_ci_median_diff(xe,xn,300)
  })
  ci <- do.call(rbind,ci)
  st[,`:=`(ci_lo=ci[,1], ci_hi=ci[,2], dataset=ds)]
}

stats <- rbindlist(lapply(levels(scores_all$dataset), stats_one_ds))
stats[, module_lab := factor(module_lab, levels=pretty_lab[MODULES])]
stats[, dataset := factor(dataset, levels=levels(scores_all$dataset))]

fwrite(stats, file.path(OUT_DIR,"ImmuneModule_stats_delta_padj_CI.tsv"), sep="\t")

L <- max(abs(stats$ci_lo), abs(stats$ci_hi), na.rm=TRUE)

plot_effect <- function(ds) {
  st <- stats[dataset==ds]
  st[, sign_col := ifelse(delta>=0, COL_UPWARM, COL_DOWNTR)]
  
  ggplot(st, aes(delta, module_lab)) +
    geom_vline(xintercept=0, linetype="dashed", color="grey40") +
    geom_errorbarh(aes(xmin=ci_lo, xmax=ci_hi), height=0.2, linewidth=0.7) +
    geom_point(aes(color=sign_col), size=4) +
    scale_color_identity() +
    geom_text(aes(label=star), nudge_x=0.09, size=7) +
    coord_cartesian(xlim=c(-L,L)) +
    theme_bw(base_size=14) +
    theme(
      axis.text.y = element_text(size=13),
      axis.title.x = element_text(size=14),
      plot.title = element_text(size=18, face="bold"),
      plot.margin = margin(10,10,10,10)
    ) +
    labs(title=ds, x="Δ immune module score (HGSC − Normal)")
}

p1 <- plot_effect(DSETS[1])
p2 <- plot_effect(DSETS[2])
p3 <- plot_effect(DSETS[3])
pM <- plot_effect("META (pooled)")

fig_effect <- (p1 | p2 | p3) / pM +
  plot_layout(heights=c(1,1.2)) +
  plot_annotation(
    title="Immune systems biology: module-level shifts in HGSC vs Normal",
    theme=theme(plot.title=element_text(size=18,face="bold"))
  )

ggsave(file.path(OUT_DIR,paste0("ImmuneModules_Effect_Lollipop_",tag,".pdf")),
       fig_effect, width=16, height=12, device=cairo_pdf)
ggsave(file.path(OUT_DIR,paste0("ImmuneModules_Effect_Lollipop_",tag,".png")),
       fig_effect, width=16, height=12, dpi=600, bg="white")

hm <- stats[,.(dataset,module,module_lab,delta,star)]
hm[, module_lab := factor(module_lab, levels=rev(pretty_lab[MODULES]))]

Lh <- max(abs(hm$delta), na.rm=TRUE)

p_hm <- ggplot(hm, aes(dataset,module_lab,fill=delta)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=paste0(sprintf("%.2f",delta),"\n",star)),
            size=5.5) +
  scale_fill_gradient2(
    low=COL_DOWNTR, mid="white", high=COL_UPWARM,
    midpoint=0, limits=c(-Lh,Lh)
  ) +
  theme_bw(base_size=14) +
  theme(
    axis.text.x = element_text(size=13, face="bold"),
    axis.text.y = element_text(size=13),
    plot.title = element_text(size=18, face="bold"),
    plot.margin = margin(10,10,10,10)
  ) +
  labs(title="Immune module delta heatmap (HGSC − Normal)")

ggsave(file.path(OUT_DIR,paste0("ImmuneModules_DeltaHeatmap_",tag,".pdf")),
       p_hm, width=11, height=8, device=cairo_pdf)
ggsave(file.path(OUT_DIR,paste0("ImmuneModules_DeltaHeatmap_",tag,".png")),
       p_hm, width=11, height=8, dpi=600, bg="white")

message("DONE. Outputs saved to: ", OUT_DIR)
