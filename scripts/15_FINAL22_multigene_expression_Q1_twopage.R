# scripts/15_FINAL22_multigene_expression_Q1_twopage.R
# Two-page PDF:
#   Page 1: Up-regulated FINAL genes
#   Page 2: Down-regulated FINAL genes
# Layout:
#   Top row: 3 datasets
#   Bottom row: META pooled

source("scripts/01_params.R")

pkgs <- c("data.table","ggplot2","patchwork","stringr")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

options(scipen=999)

DIR_MOD <- if (exists("DIR_MOD")) DIR_MOD else file.path("results","models")
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
N_GENES <- 22

COL_NORMAL <- "#9E9E9E"
COL_UP     <- "#E64B35"
COL_DOWN   <- "#00BFC4"

p_to_star <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

cand_fin_f <- file.path(DIR_MOD, paste0("STRICT_candidates_FINAL_", tag, ".tsv"))
must_exist(cand_fin_f)

final_list <- fread(cand_fin_f)

stopifnot("gene" %in% names(final_list))

if (!("meta_logFC" %in% names(final_list))) final_list[, meta_logFC := NA_real_]
if (!("meta_fdr" %in% names(final_list))) final_list[, meta_fdr := NA_real_]
if (!("n_models" %in% names(final_list))) final_list[, n_models := NA_integer_]

final_list[, direction := fifelse(is.na(meta_logFC),"NA",
                                  fifelse(meta_logFC > 0,"Up","Down"))]

setorder(final_list, meta_fdr, -n_models)

final22 <- final_list[direction %in% c("Up","Down")][1:min(N_GENES,.N)]

final_up   <- final22[direction=="Up"][order(meta_fdr)]
final_down <- final22[direction=="Down"][order(meta_fdr)]

load_expr_and_sheet <- function(ds){
  
  expr_path  <- file.path(DIR_PROC,paste0(ds,"_expr_gene.rds"))
  sheet_path <- file.path(DIR_META,paste0(ds,"_sample_sheet.csv"))
  
  must_exist(expr_path)
  must_exist(sheet_path)
  
  expr  <- readRDS(expr_path)
  sheet <- fread(sheet_path)
  
  sheet[final_group=="", final_group:=NA_character_]
  sheet <- sheet[!is.na(final_group) & final_group %in% c("HGSC","Normal")]
  
  sids <- intersect(colnames(expr),sheet$sample_id)
  
  if(length(sids) < 6) stop("Too few aligned samples for dataset: ",ds)
  
  sheet <- sheet[match(sids,sample_id)]
  expr  <- expr[,sids,drop=FALSE]
  
  list(expr=expr,sheet=sheet)
}

build_long_ds <- function(ds,genes_keep){
  
  obj <- load_expr_and_sheet(ds)
  
  expr  <- obj$expr
  sheet <- obj$sheet
  
  keep_genes <- intersect(genes_keep,rownames(expr))
  
  if(length(keep_genes) < 3) stop("Too few genes found in ",ds)
  
  mat <- expr[keep_genes,,drop=FALSE]
  
  dt <- as.data.table(t(mat),keep.rownames="sample_id")
  
  dt <- merge(dt,sheet[,.(sample_id,final_group)],by="sample_id",all.x=TRUE)
  
  setnames(dt,"final_group","group")
  
  dt <- dt[!is.na(group)]
  
  long <- melt(
    dt,
    id.vars=c("sample_id","group"),
    variable.name="gene",
    value.name="expr"
  )
  
  long[,dataset:=ds]
  long[,expr:=as.numeric(expr)]
  
  long[,expr_z:={
    m <- mean(expr,na.rm=TRUE)
    s <- sd(expr,na.rm=TRUE)
    if(!is.finite(s) || s==0) 0 else (expr-m)/s
  },by=gene]
  
  long
}

get_stats <- function(dt_one){
  
  st <- dt_one[,.(
    
    p_raw = tryCatch(
      wilcox.test(expr_z ~ group)$p.value,
      error=function(e) NA_real_
    ),
    
    y_max = max(expr_z,na.rm=TRUE)
    
  ),by=gene]
  
  st[,p_adj := p.adjust(p_raw,"BH")]
  st[,star := vapply(p_adj,p_to_star,character(1))]
  
  st[,y_pos := y_max + 0.60]
  
  st
}

plot_panel <- function(long_dt,stats_dt,ds_title,genes_levels){
  
  dd <- copy(long_dt)
  ss <- copy(stats_dt)
  
  dd[,gene := factor(gene,levels=genes_levels)]
  dd[,group := factor(group,levels=c("Normal","HGSC"))]
  
  dodge <- position_dodge(width=0.78)
  jd <- position_jitterdodge(jitter.width=0.14,dodge.width=0.78)
  
  p <- ggplot(dd,aes(x=gene,y=expr_z))
  
  p <- p +
    
    geom_violin(
      aes(fill=fill_col),
      position=dodge,
      trim=FALSE,
      alpha=0.28,
      linewidth=0.35,
      color="grey25"
    ) +
    
    geom_boxplot(
      aes(fill=fill_col),
      position=dodge,
      width=0.24,
      outlier.shape=NA,
      linewidth=0.5,
      color="black"
    ) +
    
    geom_point(
      position=jd,
      size=0.6,
      alpha=0.75,
      color="black"
    ) +
    
    geom_text(
      data=ss,
      aes(x=gene,y=y_pos,label=star),
      inherit.aes=FALSE,
      size=6,
      fontface="bold",
      vjust=0
    ) +
    
    scale_fill_identity() +
    
    theme_bw(base_size=14,base_family="Helvetica") +
    
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth=0.25,color="grey90"),
      
      plot.title = element_text(
        size=14,
        face="bold",
        margin=margin(b=6)
      ),
      
      axis.title.x = element_text(size=13,face="bold"),
      axis.title.y = element_text(size=13,face="bold"),
      
      axis.text.x = element_text(
        size=11,
        angle=45,
        hjust=1,
        vjust=1,
        face="bold"
      ),
      
      axis.text.y = element_text(size=11),
      
      plot.margin = margin(10,10,10,10)
    ) +
    
    labs(
      title=ds_title,
      x="FINAL candidate genes",
      y="Relative expression (within-dataset Z-score)"
    )
  
  p
}

make_page <- function(final_sub,page_label=c("UP","DOWN")){
  
  page_label <- match.arg(page_label)
  
  genes_keep <- final_sub$gene
  
  long_list <- lapply(DSETS,function(ds) build_long_ds(ds,genes_keep))
  names(long_list) <- DSETS
  
  long_all <- rbindlist(long_list,use.names=TRUE,fill=TRUE)
  
  dir_map <- final_sub[,.(gene,direction,meta_logFC,meta_fdr,n_models)]
  
  long_all <- merge(long_all,dir_map,by="gene",all.x=TRUE)
  
  meta_long <- copy(long_all)
  meta_long[,dataset:="META (pooled Z-score)"]
  
  long_all2 <- rbind(long_all,meta_long,use.names=TRUE,fill=TRUE)
  
  long_all2[,fill_col :=
              ifelse(group=="Normal",COL_NORMAL,
                     ifelse(direction=="Up",COL_UP,COL_DOWN))
  ]
  
  stats_list <- lapply(unique(long_all2$dataset),function(ds){
    
    dd <- long_all2[dataset==ds]
    
    s <- get_stats(dd)
    s[,dataset:=ds]
    
    s
  })
  
  stats_all <- rbindlist(stats_list,use.names=TRUE,fill=TRUE)
  
  genes_levels <- genes_keep
  
  p1 <- plot_panel(
    long_all2[dataset==DSETS[1]],
    stats_all[dataset==DSETS[1]],
    DSETS[1],
    genes_levels
  )
  
  p2 <- plot_panel(
    long_all2[dataset==DSETS[2]],
    stats_all[dataset==DSETS[2]],
    DSETS[2],
    genes_levels
  )
  
  p3 <- plot_panel(
    long_all2[dataset==DSETS[3]],
    stats_all[dataset==DSETS[3]],
    DSETS[3],
    genes_levels
  )
  
  pM <- plot_panel(
    long_all2[dataset=="META (pooled Z-score)"],
    stats_all[dataset=="META (pooled Z-score)"],
    "META (pooled Z-score)",
    genes_levels
  )
  
  fig <- (p1 | p2 | p3) / pM +
    plot_layout(heights=c(1.0,1.35)) +
    plot_annotation(
      title=sprintf("FINAL candidates (%s) | n=%d",page_label,length(genes_levels)),
      theme=theme(
        plot.title=element_text(size=16,face="bold",family="Helvetica")
      )
    )
  
  list(fig=fig,genes=genes_levels)
}

page_up   <- make_page(final_up,"UP")
page_down <- make_page(final_down,"DOWN")

OUT_PDF <- file.path(
  DIR_FIG,
  paste0("FINAL22_TwoPage_UpDown_BoxViolin_Meta_",tag,".pdf")
)

OUT_PNG_UP <- file.path(
  DIR_FIG,
  paste0("FINAL22_UP_BoxViolin_Meta_",tag,".png")
)

OUT_PNG_DN <- file.path(
  DIR_FIG,
  paste0("FINAL22_DOWN_BoxViolin_Meta_",tag,".png")
)

grDevices::cairo_pdf(OUT_PDF,width=16.5,height=11.8)
print(page_up$fig)
print(page_down$fig)
dev.off()

ggsave(
  OUT_PNG_UP,
  page_up$fig,
  width=16.5,
  height=11.8,
  units="in",
  dpi=600,
  bg="white"
)

ggsave(
  OUT_PNG_DN,
  page_down$fig,
  width=16.5,
  height=11.8,
  units="in",
  dpi=600,
  bg="white"
)

fwrite(
  final_up,
  file.path(DIR_FIG,paste0("FINAL_UP_gene_list_",tag,".tsv")),
  sep="\t"
)

fwrite(
  final_down,
  file.path(DIR_FIG,paste0("FINAL_DOWN_gene_list_",tag,".tsv")),
  sep="\t"
)

message("Saved two-page PDF: ",OUT_PDF)
