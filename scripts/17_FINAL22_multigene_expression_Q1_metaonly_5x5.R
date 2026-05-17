# FINAL 22 genes | META pooled only
# Journal-grade (Q1)
# Violin + Box + Enhanced sample points + BH-adjusted Wilcoxon stars

source("scripts/01_params.R")

pkgs <- c("data.table","ggplot2","stringr","cowplot")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(cowplot)
})

options(scipen=999)

DIR_MOD  <- if (exists("DIR_MOD"))  DIR_MOD  else file.path("results","models")
DIR_FIG  <- if (exists("DIR_FIG"))  DIR_FIG  else "figures"
DIR_PROC <- if (exists("DIR_PROC")) DIR_PROC else file.path("results","processed")
DIR_META <- if (exists("DIR_META")) DIR_META else "metadata"

dir.create(DIR_FIG, showWarnings=FALSE, recursive=TRUE)

must_exist <- function(f) if (!file.exists(f)) stop("Missing file: ", f)

safe_tag <- function(x) {
  x <- gsub("<=|≥|≤", "le", x)
  x <- gsub("<", "lt", x)
  x <- gsub(">", "gt", x)
  x <- gsub("=", "eq", x)
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

# ---------------- SETTINGS ----------------
TOPK <- 50
MIN_MODELS <- 3
META_FDR_CUT <- 0.05

tag <- safe_tag(sprintf(
  "TOPK%d_minModels%d_metaFDRle%.2f",
  TOPK, MIN_MODELS, META_FDR_CUT
))

DSETS <- c("GSE14407","GSE38666","GSE52037")

N_GENES <- 22
N_COLS  <- 5

COL_NORMAL <- "#9E9E9E"
COL_UP     <- "#D84A3A"
COL_DOWN   <- "#14B3C7"
STRIP_FILL <- "#D9D9D9"

p_to_star <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# ---------------- INPUT ----------------
cand_fin_f <- file.path(
  DIR_MOD,
  paste0("STRICT_candidates_FINAL_", tag, ".tsv")
)
must_exist(cand_fin_f)

final_list <- fread(cand_fin_f)
final_list[, direction := fifelse(meta_logFC > 0, "Up", "Down")]

setorder(final_list, meta_fdr, -n_models)

final22 <- final_list[
  direction %in% c("Up","Down")
][1:min(N_GENES, .N)]

final22[, direction := factor(direction, levels=c("Up","Down"))]
setorder(final22, direction, meta_fdr, -n_models)

genes_levels <- final22$gene

n_up <- sum(final22$direction=="Up", na.rm=TRUE)
n_dn <- sum(final22$direction=="Down", na.rm=TRUE)

# ---------------- DATA LOAD ----------------
load_expr_and_sheet <- function(ds) {
  
  expr  <- readRDS(file.path(DIR_PROC, paste0(ds, "_expr_gene.rds")))
  sheet <- fread(file.path(DIR_META, paste0(ds, "_sample_sheet.csv")))
  
  sheet <- sheet[final_group %in% c("HGSC","Normal")]
  
  sids <- intersect(colnames(expr), sheet$sample_id)
  sheet <- sheet[match(sids, sample_id)]
  expr  <- expr[, sids, drop=FALSE]
  
  list(expr=expr, sheet=sheet)
}

build_long_ds <- function(ds, genes_keep) {
  
  obj <- load_expr_and_sheet(ds)
  expr  <- obj$expr
  sheet <- obj$sheet
  
  mat <- expr[intersect(genes_keep, rownames(expr)), , drop=FALSE]
  
  dt <- as.data.table(t(mat), keep.rownames="sample_id")
  dt <- merge(dt, sheet[,.(sample_id,final_group)], by="sample_id")
  setnames(dt,"final_group","group")
  
  long <- melt(
    dt,
    id.vars=c("sample_id","group"),
    variable.name="gene",
    value.name="expr"
  )
  
  long[, dataset:=ds]
  
  long[, expr_z := (expr - mean(expr,na.rm=TRUE)) /
         sd(expr,na.rm=TRUE),
       by=gene]
  
  long
}

# ---------------- META MERGE ----------------
meta_long <- rbindlist(
  lapply(DSETS, build_long_ds, genes_keep=genes_levels)
)

meta_long <- merge(
  meta_long,
  final22[,.(gene,direction)],
  by="gene"
)

meta_long[, group:=factor(group, levels=c("Normal","HGSC"))]
meta_long[, gene :=factor(gene,  levels=genes_levels)]

meta_long[, fill_col := ifelse(
  group=="Normal",
  COL_NORMAL,
  ifelse(direction=="Up", COL_UP, COL_DOWN)
)]

# ---------------- STATS ----------------
st <- meta_long[, .(
  p_raw = tryCatch(
    wilcox.test(expr_z ~ group)$p.value,
    error=function(e) NA_real_
  ),
  y_max = max(expr_z, na.rm=TRUE)
), by=gene]

st[, p_adj := p.adjust(p_raw, "BH")]
st[, star  := vapply(p_adj, p_to_star, character(1))]
st[, y_pos := y_max + 0.6]

# ---------------- PLOT ----------------
jd <- position_jitter(width=0.08, height=0, seed=123)

p_core <- ggplot(meta_long, aes(group, expr_z)) +
  
  geom_violin(
    aes(fill=fill_col),
    trim=FALSE,
    alpha=0.18,
    linewidth=0.35,
    color="grey30"
  ) +
  
  geom_boxplot(
    aes(fill=fill_col),
    width=0.48,
    outlier.shape=NA,
    linewidth=1.05,
    color="black",
    fatten=3.8,
    alpha=0.80
  ) +
  
  geom_point(
    aes(fill=fill_col),
    position=jd,
    shape=21,
    size=1.55,
    stroke=0.35,
    color="black",
    alpha=0.95
  ) +
  
  geom_text(
    data=st,
    aes(x="HGSC", y=y_pos, label=star),
    inherit.aes=FALSE,
    size=5.5,
    fontface="bold"
  ) +
  
  scale_fill_identity() +
  
  facet_wrap(~gene, ncol=N_COLS) +
  
  theme_bw(base_size=13) +
  
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="#EAEAEA",linewidth=0.3),
    
    strip.text=element_text(size=11, face="italic"),
    strip.background=element_rect(
      fill=STRIP_FILL,
      color="grey60",
      linewidth=0.35
    ),
    
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=12),
    axis.text=element_text(size=10),
    
    axis.line=element_line(linewidth=0.35),
    axis.ticks=element_line(linewidth=0.35),
    
    plot.title=element_text(
      size=16,
      face="bold",
      margin=ggplot2::margin(b=6, unit="pt")
    ),
    
    plot.subtitle=element_text(
      size=11,
      margin=ggplot2::margin(b=4, unit="pt")
    ),
    
    plot.margin=ggplot2::margin(14,14,6,14,unit="pt")
  ) +
  
  labs(
    y="Relative expression (within-dataset Z-score; pooled META)",
    title=sprintf(
      "FINAL 22 candidate genes (META pooled) | Up=%d, Down=%d",
      n_up, n_dn
    ),
    subtitle="Points = individual samples | Stars = BH-adjusted Wilcoxon"
  )

y_bottom <- min(meta_long$expr_z, na.rm=TRUE)
y_top    <- max(st$y_pos, na.rm=TRUE)

p_core <- p_core +
  coord_cartesian(
    ylim=c(y_bottom, y_top + 0.4),
    clip="off"
  )

# ---------------- EXPORT (3 HEIGHT MODES) ----------------

# A4
OUT_A4_PDF <- file.path(DIR_FIG,
                        paste0("FINAL22_METAonly_5x5_", tag, "_A4.pdf"))
OUT_A4_PNG <- file.path(DIR_FIG,
                        paste0("FINAL22_METAonly_5x5_", tag, "_A4.png"))

cairo_pdf(OUT_A4_PDF, width=12, height=12)
print(p_core)
dev.off()

ggsave(OUT_A4_PNG, p_core,
       width=12, height=12,
       units="in", dpi=700,
       bg="white", limitsize=FALSE)

# TALL
OUT_TALL_PDF <- file.path(DIR_FIG,
                          paste0("FINAL22_METAonly_5x5_", tag, "_TALL.pdf"))
OUT_TALL_PNG <- file.path(DIR_FIG,
                          paste0("FINAL22_METAonly_5x5_", tag, "_TALL.png"))

cairo_pdf(OUT_TALL_PDF, width=12, height=16)
print(p_core)
dev.off()

ggsave(OUT_TALL_PNG, p_core,
       width=12, height=16,
       units="in", dpi=700,
       bg="white", limitsize=FALSE)

# XL
OUT_XL_PDF <- file.path(DIR_FIG,
                        paste0("FINAL22_METAonly_5x5_", tag, "_XL.pdf"))
OUT_XL_PNG <- file.path(DIR_FIG,
                        paste0("FINAL22_METAonly_5x5_", tag, "_XL.png"))

cairo_pdf(OUT_XL_PDF, width=12, height=20)
print(p_core)
dev.off()

ggsave(OUT_XL_PNG, p_core,
       width=12, height=20,
       units="in", dpi=700,
       bg="white", limitsize=FALSE)

# Save gene list
fwrite(final22,
       file.path(DIR_FIG,
                 paste0("FINAL22_gene_list_", tag, ".tsv")),
       sep="\t"
)

message("✅ FINAL22 figure exported successfully.")
