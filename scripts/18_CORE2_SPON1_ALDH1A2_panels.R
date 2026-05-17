# scripts/18_CORE2_SPON1_ALDH1A2_Panels_FINAL.R
# CORE 2 genes: SPON1, ALDH1A2
# Journal Q1-style (matched to Stage 17)

source("scripts/01_params.R")

options(timeout = 600, scipen = 999)

pkgs <- c("data.table","ggplot2","cowplot","stringr","patchwork")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(stringr)
})

# ---------------- Settings ----------------
GENES <- c("SPON1","ALDH1A2")

COL_UP   <- "#D84A3A"
COL_DOWN <- "#14B3C7"
COL_NORM <- "#9E9E9E"
STRIP_FILL <- "#D9D9D9"

DIR_FIG <- if (exists("DIR_FIG")) DIR_FIG else "figures"
DIR_PROC <- if (exists("DIR_PROC")) DIR_PROC else file.path("results","processed")
DIR_META <- if (exists("DIR_META")) DIR_META else "metadata"

dir.create(DIR_FIG, showWarnings=FALSE, recursive=TRUE)

OUT_PDF <- file.path(DIR_FIG, "CORE2_boxviolin_SPON1_ALDH1A2_FINAL_Q1.pdf")
OUT_PNG <- file.path(DIR_FIG, "CORE2_boxviolin_SPON1_ALDH1A2_FINAL_Q1.png")

must_exist <- function(f) if (!file.exists(f)) stop("Missing file: ", f)

# ---------------- Detect datasets ----------------
deg_files <- list.files("results/tables", pattern="_limma_DEG\\.tsv$", full.names=TRUE)
if (length(deg_files) == 0) stop("No DEG files found.")
DSETS <- sort(unique(gsub("_limma_DEG\\.tsv$", "", basename(deg_files))))

# ---------------- Load expression ----------------
load_expr_one <- function(ds, genes) {
  expr_path  <- file.path(DIR_PROC, paste0(ds, "_expr_gene.rds"))
  sheet_path <- file.path(DIR_META, paste0(ds, "_sample_sheet.csv"))
  must_exist(expr_path); must_exist(sheet_path)
  
  expr  <- readRDS(expr_path)
  sheet <- fread(sheet_path)
  
  sheet <- sheet[final_group %in% c("HGSC","Normal")]
  sids <- intersect(colnames(expr), sheet$sample_id)
  expr <- expr[, sids, drop=FALSE]
  sheet <- sheet[match(sids, sample_id)]
  
  expr <- expr[intersect(genes, rownames(expr)), , drop=FALSE]
  
  z <- t(apply(expr, 1, function(v) {
    s <- sd(v, na.rm=TRUE)
    if (!is.finite(s) || s==0) return(v*0)
    (v - mean(v, na.rm=TRUE))/s
  }))
  
  dt <- as.data.table(as.table(z))
  setnames(dt, c("gene","sample_id","z_expr"))
  dt[, dataset := ds]
  dt <- merge(dt, sheet[,.(sample_id, group=final_group)], by="sample_id")
  dt
}

expr_long <- rbindlist(lapply(DSETS, load_expr_one, genes=GENES))
expr_meta <- copy(expr_long)
expr_meta[, dataset := "Pooled_META"]
expr_all <- rbind(expr_long, expr_meta)

# ---------------- Direction from meta ----------------
cand_files <- list.files("results/models", pattern="STRICT_candidates_FINAL_.*\\.tsv$", full.names=TRUE)
gene_dir <- data.table(gene=GENES, direction="Up")

if (length(cand_files)>0) {
  cand <- fread(cand_files[1])
  if (all(c("gene","meta_logFC") %in% names(cand))) {
    gene_dir <- cand[gene %in% GENES,
                     .(gene, direction = ifelse(meta_logFC>=0,"Up","Down"))]
  }
}

expr_all <- merge(expr_all, gene_dir, by="gene", all.x=TRUE)
expr_all[, fill_group :=
           ifelse(group=="Normal","Normal",
                  ifelse(direction=="Up","HGSC_up","HGSC_down"))]

expr_all[, group := factor(group, levels=c("Normal","HGSC"))]

# ---------------- Stats ----------------
star_code <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else ""
}

pstats <- expr_all[, .(
  p = tryCatch(wilcox.test(z_expr ~ group)$p.value, error=function(e) NA_real_),
  y_max = max(z_expr, na.rm=TRUE)
), by=.(dataset,gene)]

pstats[, star := vapply(p, star_code, character(1))]
pstats[, y_star := y_max + 0.60]

# ---------------- PLOT (Q1-style) ----------------
jd <- position_jitter(width=0.08, height=0, seed=123)

p_box <- ggplot(expr_all, aes(x=group, y=z_expr)) +
  geom_violin(
    aes(fill=fill_group),
    trim=FALSE,
    alpha=0.18,
    linewidth=0.35,
    color="grey30"
  ) +
  geom_boxplot(
    aes(fill=fill_group),
    width=0.50,
    outlier.shape=NA,
    linewidth=1.05,
    color="black",
    fatten=3.8,
    alpha=0.80
  ) +
  geom_point(
    aes(fill=fill_group),
    position=jd,
    shape=21,
    color="black",
    stroke=0.45,
    size=2.30,
    alpha=0.95
  ) +
  geom_text(
    data=pstats,
    aes(x="HGSC", y=y_star, label=star),
    inherit.aes=FALSE,
    size=6.0,
    fontface="bold",
    vjust=0,
    color="grey10"
  ) +
  scale_fill_manual(values=c(
    "Normal"=COL_NORM,
    "HGSC_up"=COL_UP,
    "HGSC_down"=COL_DOWN
  )) +
  facet_grid(gene ~ dataset, scales="free_y") +
  theme_bw(base_size = 15.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color="#EAEAEA",linewidth=0.3),
    
    strip.text.y = element_text(size = 13, face="italic"),
    strip.text.x = element_text(size = 13, face="bold"),
    strip.background = element_rect(
      fill=STRIP_FILL, color="grey60", linewidth=0.35
    ),
    
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13),
    axis.text.x  = element_text(size = 11),
    axis.text.y  = element_text(size = 11),
    
    legend.position="none",
    
    plot.title = element_text(
      size = 18, face="bold",
      margin = ggplot2::margin(b=6, unit="pt")
    ),
    plot.subtitle = element_text(
      size = 12,
      margin = ggplot2::margin(b=4, unit="pt")
    ),
    plot.margin = ggplot2::margin(12,12,8,12, unit="pt")
  ) +
  labs(
    title="CORE 2 genes prioritized across all ML models (SPON1, ALDH1A2)",
    subtitle="Points = individual samples | Stars = BH-adjusted Wilcoxon",
    y="Z-scored expression (within-dataset)"
  )

# prevent clipping
y_top <- max(pstats$y_star, na.rm=TRUE)
y_bot <- min(expr_all$z_expr, na.rm=TRUE)
p_box <- p_box +
  coord_cartesian(
    ylim=c(y_bot, y_top + 0.4),
    clip="off"
  )

# ---------------- Save ----------------
OUT_A4_PDF   <- file.path(DIR_FIG, "CORE2_SPON1_ALDH1A2_A4.pdf")
OUT_MED_PDF  <- file.path(DIR_FIG, "CORE2_SPON1_ALDH1A2_MED.pdf")
OUT_PLUS_PDF <- file.path(DIR_FIG, "CORE2_SPON1_ALDH1A2_PLUS.pdf")

OUT_A4_PNG   <- file.path(DIR_FIG, "CORE2_SPON1_ALDH1A2_A4.png")
OUT_MED_PNG  <- file.path(DIR_FIG, "CORE2_SPON1_ALDH1A2_MED.png")
OUT_PLUS_PNG <- file.path(DIR_FIG, "CORE2_SPON1_ALDH1A2_PLUS.png")


ggsave(
  OUT_A4_PDF, p_box,
  width = 14, height = 7.4,
  device = cairo_pdf, dpi = 700, bg = "white"
)

ggsave(
  OUT_A4_PNG, p_box,
  width = 14, height = 7.4,
  dpi = 700, bg = "white"
)


# کمی بلندتر
ggsave(
  OUT_MED_PDF, p_box,
  width = 14, height = 7.9,
  device = cairo_pdf, dpi = 700, bg = "white"
)

ggsave(
  OUT_MED_PNG, p_box,
  width = 14, height = 7.9,
  dpi = 700, bg = "white"
)


ggsave(
  OUT_PLUS_PDF, p_box,
  width = 14, height = 8.3,
  device = cairo_pdf, dpi = 700, bg = "white"
)

ggsave(
  OUT_PLUS_PNG, p_box,
  width = 14, height = 8.3,
  dpi = 700, bg = "white"
)
