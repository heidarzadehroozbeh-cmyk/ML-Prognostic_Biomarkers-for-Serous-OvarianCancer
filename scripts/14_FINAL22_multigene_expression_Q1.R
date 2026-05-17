# scripts/14_FINAL22_multigene_expression_Q1.R
# Single-page figure: 22 FINAL genes (Up/Down) with Box+Violin+Points, significance stars
# Layout: 3 datasets on first row (3 columns) + META pooled on second row (full width)
# Colors: Normal=grey; HGSC-Up=warm; HGSC-Down=turquoise

source("scripts/01_params.R")

pkgs <- c("data.table","ggplot2","patchwork","stringr")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(stringr)
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
if (!("meta_fdr"   %in% names(final_list))) final_list[, meta_fdr   := NA_real_]
if (!("n_models"   %in% names(final_list))) final_list[, n_models   := NA_integer_]

final_list[, direction := fifelse(is.na(meta_logFC), "NA",
                                  fifelse(meta_logFC > 0, "Up", "Down"))]

setorder(final_list, meta_fdr, -n_models)
final22 <- final_list[direction %in% c("Up","Down")][1:min(N_GENES, .N)]

final22 <- rbind(
  final22[direction=="Up"][order(meta_fdr)],
  final22[direction=="Down"][order(meta_fdr)]
)

genes <- final22$gene

load_expr_and_sheet <- function(ds) {
  
  expr_path  <- file.path(DIR_PROC, paste0(ds, "_expr_gene.rds"))
  sheet_path <- file.path(DIR_META, paste0(ds, "_sample_sheet.csv"))
  
  must_exist(expr_path)
  must_exist(sheet_path)
  
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

build_long_ds <- function(ds) {
  
  obj <- load_expr_and_sheet(ds)
  expr <- obj$expr
  sheet <- obj$sheet
  
  keep_genes <- intersect(genes, rownames(expr))
  if (length(keep_genes) < 5) stop("Too few selected genes found in ", ds)
  
  mat <- expr[keep_genes, , drop=FALSE]
  
  dt <- as.data.table(t(mat), keep.rownames="sample_id")
  
  dt <- merge(dt, sheet[, .(sample_id, final_group)], by="sample_id", all.x=TRUE)
  setnames(dt, "final_group", "group")
  dt <- dt[!is.na(group)]
  
  long <- melt(dt,
               id.vars=c("sample_id","group"),
               variable.name="gene",
               value.name="expr")
  
  long[, dataset := ds]
  
  long[, expr := as.numeric(expr)]
  
  long[, expr_z := {
    m <- mean(expr, na.rm=TRUE)
    s <- sd(expr, na.rm=TRUE)
    if (!is.finite(s) || s==0) 0 else (expr - m)/s
  }, by=gene]
  
  long
}

long_list <- lapply(DSETS, build_long_ds)
names(long_list) <- DSETS
long_all <- rbindlist(long_list, use.names=TRUE, fill=TRUE)

dir_map <- final22[, .(gene, direction, meta_logFC, meta_fdr, n_models)]
long_all <- merge(long_all, dir_map, by="gene", all.x=TRUE)

meta_long <- copy(long_all)
meta_long[, dataset := "META (pooled Z-score)"]

long_all2 <- rbind(long_all, meta_long, use.names=TRUE, fill=TRUE)

long_all2[, gene := factor(gene, levels=genes)]
long_all2[, group := factor(group, levels=c("Normal","HGSC"))]

long_all2[, fill_col :=
            ifelse(group=="Normal",
                   COL_NORMAL,
                   ifelse(direction=="Up", COL_UP, COL_DOWN))]

get_stats <- function(dt_one) {
  
  st <- dt_one[, .(
    p_raw = tryCatch(wilcox.test(expr_z ~ group)$p.value,
                     error=function(e) NA_real_),
    y_max = max(expr_z, na.rm=TRUE)
  ), by=gene]
  
  st[, p_adj := p.adjust(p_raw, method="BH")]
  st[, star := vapply(p_adj, p_to_star, character(1))]
  st[, y_pos := y_max + 0.35]
  
  st
}

stats_list <- lapply(unique(long_all2$dataset), function(ds) {
  
  dd <- long_all2[dataset==ds]
  
  if (nrow(dd)==0) return(NULL)
  
  s <- get_stats(dd)
  s[, dataset := ds]
  
  s
})

stats_all <- rbindlist(stats_list, use.names=TRUE, fill=TRUE)

plot_panel <- function(ds, show_x=TRUE) {
  
  dd <- long_all2[dataset==ds]
  ss <- stats_all[dataset==ds]
  
  dodge <- position_dodge(width=0.75)
  jd <- position_jitterdodge(jitter.width=0.12, dodge.width=0.75)
  
  p <- ggplot(dd, aes(x=gene, y=expr_z)) +
    
    geom_violin(aes(fill=fill_col),
                position=dodge,
                trim=FALSE,
                alpha=0.22,
                linewidth=0.25,
                color="grey20") +
    
    geom_boxplot(aes(fill=fill_col),
                 position=dodge,
                 width=0.22,
                 outlier.shape=NA,
                 linewidth=0.35,
                 color="grey15") +
    
    geom_point(position=jd,
               size=0.35,
               alpha=0.75,
               color="black") +
    
    geom_text(data=ss,
              aes(x=gene, y=y_pos, label=star),
              inherit.aes=FALSE,
              size=3.2,
              vjust=0) +
    
    scale_fill_identity() +
    
    theme_bw(base_size=11) +
    
    theme(
      panel.grid.minor = element_blank(),
      
      plot.title = element_text(
        size=12,
        face="bold",
        margin=ggplot2::margin(b=6, unit="pt")
      ),
      
      axis.title.x = element_text(size=11),
      axis.title.y = element_text(size=11),
      
      axis.text.x = element_text(
        size=7,
        angle=45,
        hjust=1,
        vjust=1
      ),
      
      axis.text.y = element_text(size=9),
      
      plot.margin = ggplot2::margin(8,8,8,8,"pt")
    ) +
    
    labs(
      title = ds,
      x = if (show_x)
        "FINAL candidate genes (Up/Down ordered by metaFDR)"
      else NULL,
      y = "Relative expression (within-dataset Z-score)"
    )
  
  if (!show_x) {
    
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
    
  }
  
  p
}

p1 <- plot_panel(DSETS[1], show_x=FALSE)
p2 <- plot_panel(DSETS[2], show_x=FALSE)
p3 <- plot_panel(DSETS[3], show_x=FALSE)

pM <- plot_panel("META (pooled Z-score)", show_x=TRUE) +
  
  plot_annotation(
    subtitle =
      "Normal = grey; HGSC-Up = warm; HGSC-Down = turquoise. Points = individual samples. Stars = BH-adjusted Wilcoxon (* <0.05, ** <0.01, *** <0.001).",
    theme = theme(
      plot.subtitle = element_text(size=10)
    )
  )

fig <- (p1 | p2 | p3) / pM +
  
  plot_layout(
    heights=c(1.0, 1.25)
  ) +
  
  plot_annotation(
    title =
      sprintf(
        "FINAL candidates (n=%d): multi-dataset relative expression + pooled meta panel",
        length(genes)
      ),
    theme = theme(
      plot.title = element_text(
        size=14,
        face="bold"
      )
    )
  )

OUT_PDF <- file.path(
  DIR_FIG,
  paste0(
    "FINAL22_MultiGene_BoxViolin_Meta_",
    tag,
    ".pdf"
  )
)

OUT_PNG <- file.path(
  DIR_FIG,
  paste0(
    "FINAL22_MultiGene_BoxViolin_Meta_",
    tag,
    ".png"
  )
)

ggsave(
  OUT_PDF,
  fig,
  width=14.5,
  height=11.5,
  units="in",
  device=cairo_pdf,
  limitsize=FALSE
)

ggsave(
  OUT_PNG,
  fig,
  width=14.5,
  height=11.5,
  units="in",
  dpi=600,
  bg="white",
  limitsize=FALSE
)

fwrite(
  final22,
  file.path(
    DIR_FIG,
    paste0(
      "FINAL22_gene_list_",
      tag,
      ".tsv"
    )
  ),
  sep="\t"
)

message("Saved: ", OUT_PDF)
