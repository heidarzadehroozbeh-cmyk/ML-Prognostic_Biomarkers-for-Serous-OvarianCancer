# scripts/07_publication_figures_pack_v3.R
# Publication figure pack (v3)
# Improvements:
# - Heatmaps exported as PDF + PNG (300 dpi)
# - Safe device wrappers
# - Hard file checks
# - Publication ready figures

source("scripts/01_params.R")

# ---- packages ----
pkgs <- c("data.table","ggplot2","pROC","glmnet","pheatmap","metafor")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pROC)
  library(glmnet)
  library(pheatmap)
  library(metafor)
})

# ---- helpers ----
safe_dir <- function(path) dir.create(path, showWarnings = FALSE, recursive = TRUE)

safe_pdf <- function(file, width=10, height=7, expr) {
  pdf(file, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

safe_png <- function(file, width=10, height=7, res=300, expr) {
  png(file, width = width, height = height, units="in", res=res)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

must_exist <- function(f) {
  ok <- file.exists(f) && isTRUE(file.info(f)$size > 0)
  if (!ok) stop("File not created or empty: ", f)
}

# ---- directories ----
safe_dir("figures")
safe_dir("figures/roc")
safe_dir("figures/heatmap")
safe_dir("figures/volcano")
safe_dir("figures/panel")
safe_dir("figures/forest")

# ---- knobs ----
TOP_PANEL_GENES_HEATMAP <- 40
TOP_META_UP <- 25
TOP_META_DOWN <- 25

# ================================
# 1) Load panel genes
# ================================
pan_file <- file.path(DIR_MOD, "LODO_nested_panels.tsv")
if (!file.exists(pan_file)) stop("Missing: ", pan_file)

pan <- fread(pan_file)

support <- pan[, .(
  support_folds = uniqueN(test_dataset),
  mean_freq = mean(freq, na.rm=TRUE),
  max_freq  = max(freq, na.rm=TRUE)
), by=gene][order(-support_folds,-mean_freq,-max_freq)]

final_majority <- support[support_folds >= 2]
final_strict <- support[support_folds == length(GSE_mrna)]

write.table(final_majority,file.path(DIR_MOD,"FINAL_panel_majority.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)

write.table(final_strict,file.path(DIR_MOD,"FINAL_panel_strict.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)

panel <- if (nrow(final_strict) >= 10) final_strict else final_majority
panel_genes_top <- head(panel$gene, TOP_PANEL_GENES_HEATMAP)

# ================================
# 2) Dataset loader
# ================================
load_dataset <- function(gse_id, genes_keep){
  
  expr_gene <- readRDS(file.path(DIR_PROC,paste0(gse_id,"_expr_gene.rds")))
  sheet <- read.csv(file.path(DIR_META,paste0(gse_id,"_sample_sheet.csv")))
  
  sheet$final_group[sheet$final_group==""] <- NA
  sheet <- sheet[!is.na(sheet$final_group) &
                   sheet$final_group %in% c("HGSC","Normal"),]
  
  genes <- intersect(genes_keep, rownames(expr_gene))
  X <- t(expr_gene[genes,sheet$sample_id,drop=FALSE])
  
  list(
    X=X,
    group=sheet$final_group,
    ds=rep(gse_id,nrow(X)),
    sample_id=sheet$sample_id
  )
  
}

zscore <- function(x) scale(x)

# ================================
# 3) Heatmap per dataset
# ================================
make_heatmap_one <- function(ds_id, genes_target){
  
  d <- load_dataset(ds_id, genes_target)
  
  genes_use <- intersect(genes_target, colnames(d$X))
  X <- d$X[,genes_use,drop=FALSE]
  
  ord <- order(d$group)
  X <- X[ord,,drop=FALSE]
  grp <- d$group[ord]
  
  Xz <- zscore(X)
  
  anno <- data.frame(Group=grp)
  rownames(anno) <- rownames(Xz)
  
  draw_plot <- function(){
    pheatmap::pheatmap(
      t(Xz),
      annotation_col=anno,
      show_colnames=FALSE,
      fontsize_row=8,
      main=paste0(ds_id," panel heatmap (",length(genes_use)," genes)")
    )
  }
  
  pdf_file <- file.path("figures/heatmap",paste0("Heatmap_panel_",ds_id,".pdf"))
  png_file <- file.path("figures/heatmap",paste0("Heatmap_panel_",ds_id,".png"))
  
  safe_pdf(pdf_file,10,7,draw_plot())
  safe_png(png_file,10,7,300,draw_plot())
  
  must_exist(pdf_file)
  must_exist(png_file)
  
}

for(ds in GSE_mrna) make_heatmap_one(ds,panel_genes_top)

# ================================
# 4) Combined heatmap
# ================================
hm_list <- lapply(GSE_mrna,function(x) load_dataset(x,panel_genes_top))

genes_common <- Reduce(intersect,lapply(hm_list,function(x) colnames(x$X)))

hm_X <- do.call(rbind,lapply(hm_list,function(x) x$X[,genes_common,drop=FALSE]))
hm_group <- unlist(lapply(hm_list,function(x) x$group))
hm_ds <- unlist(lapply(hm_list,function(x) x$ds))

ord <- order(hm_ds,hm_group)

hm_X <- hm_X[ord,,drop=FALSE]
hm_group <- hm_group[ord]
hm_ds <- hm_ds[ord]

for(d in unique(hm_ds)){
  idx <- which(hm_ds==d)
  hm_X[idx,] <- scale(hm_X[idx,,drop=FALSE])
}

anno_all <- data.frame(Dataset=hm_ds,Group=hm_group)
rownames(anno_all) <- rownames(hm_X)

draw_all <- function(){
  pheatmap::pheatmap(
    t(hm_X),
    annotation_col=anno_all,
    show_colnames=FALSE,
    fontsize_row=8,
    main=paste0("Panel heatmap ALL datasets (",length(genes_common)," genes)")
  )
}

pdf_all <- "figures/heatmap/Heatmap_panel_ALL_datasets.pdf"
png_all <- "figures/heatmap/Heatmap_panel_ALL_datasets.png"

safe_pdf(pdf_all,12,8,draw_all())
safe_png(png_all,12,8,300,draw_all())

must_exist(pdf_all)
must_exist(png_all)

# ================================
# 5) Meta heatmap
# ================================
meta_file <- file.path(DIR_TAB,"META_all_genes.tsv")

if(file.exists(meta_file)){
  
  meta <- fread(meta_file)
  
  meta_sig <- meta[order(meta$meta_FDR)]
  
  up <- meta_sig[meta_logFC>0][1:TOP_META_UP]
  dn <- meta_sig[meta_logFC<0][1:TOP_META_DOWN]
  
  meta_genes <- unique(c(up$gene,dn$gene))
  
  mh_list <- lapply(GSE_mrna,function(x) load_dataset(x,meta_genes))
  
  genes_common <- Reduce(intersect,lapply(mh_list,function(x) colnames(x$X)))
  
  mh_X <- do.call(rbind,lapply(mh_list,function(x) x$X[,genes_common,drop=FALSE]))
  mh_group <- unlist(lapply(mh_list,function(x) x$group))
  mh_ds <- unlist(lapply(mh_list,function(x) x$ds))
  
  ord <- order(mh_ds,mh_group)
  
  mh_X <- mh_X[ord,,drop=FALSE]
  mh_group <- mh_group[ord]
  mh_ds <- mh_ds[ord]
  
  for(d in unique(mh_ds)){
    idx <- which(mh_ds==d)
    mh_X[idx,] <- scale(mh_X[idx,,drop=FALSE])
  }
  
  anno <- data.frame(Dataset=mh_ds,Group=mh_group)
  rownames(anno) <- rownames(mh_X)
  
  draw_meta <- function(){
    pheatmap::pheatmap(
      t(mh_X),
      annotation_col=anno,
      show_colnames=FALSE,
      fontsize_row=8,
      main=paste0("META heatmap (",length(genes_common)," genes)")
    )
  }
  
  pdf_meta <- "figures/heatmap/Heatmap_META_topUpDown_ALL.pdf"
  png_meta <- "figures/heatmap/Heatmap_META_topUpDown_ALL.png"
  
  safe_pdf(pdf_meta,12,8,draw_meta())
  safe_png(png_meta,12,8,300,draw_meta())
  
  must_exist(pdf_meta)
  must_exist(png_meta)
  
}

# ================================
# 6) Volcano plots
# ================================
make_volcano <- function(gse){
  
  tt <- fread(file.path(DIR_TAB,paste0(gse,"_limma_all.tsv")))
  
  tt$neglog10 <- -log10(tt$P.Value)
  
  p <- ggplot(tt,aes(logFC,neglog10))+
    geom_point(size=1)+
    theme_bw()+
    labs(title=paste("Volcano",gse),
         x="log2FC",y="-log10(P)")
  
  pdf_file <- file.path("figures/volcano",paste0("Volcano_",gse,".pdf"))
  ggsave(pdf_file,p,width=6,height=5)
  
}

for(ds in GSE_mrna) make_volcano(ds)

# ================================
# 7) Manifest
# ================================
files <- list.files("figures",recursive=TRUE,full.names=TRUE)

man <- data.frame(
  file=files,
  bytes=file.info(files)$size
)

write.csv(man,"figures/_MANIFEST.csv",row.names=FALSE)

cat("\nFigures created:",nrow(man),"\n")
