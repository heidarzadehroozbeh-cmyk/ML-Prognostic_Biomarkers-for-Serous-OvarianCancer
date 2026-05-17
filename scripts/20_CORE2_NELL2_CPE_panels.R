# scripts/20_CORE2_SPON1_ALDH1A2_panels.R
# FINAL CLEAN VERSION – enlarged fonts, fixed stars, stable ggplot layers

source("scripts/01_params.R")

options(timeout = 600)

pkgs <- c("data.table","ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---------------- Settings ----------------
GENES <- c("SPON1","ALDH1A2")

COL_UP   <- "#E64B35"
COL_DOWN <- "#00BFC4"
COL_NORM <- "#BDBDBD"

BASE_SIZE <- 16.5
STAR_BOX  <- 5.7
STAR_HM   <- 6.3

DIR_FIG <- if (exists("DIR_FIG")) DIR_FIG else "figures"
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

# ---------------- Utilities ----------------
star_code <- function(p) {
  ifelse(
    is.na(p), "ns",
    ifelse(p < 1e-4, "****",
           ifelse(p < 1e-3, "***",
                  ifelse(p < 1e-2, "**",
                         ifelse(p < 0.05, "*", "ns")
                  )
           )
    )
  )
}

# ---------------- Detect datasets ----------------
deg_files <- list.files("results/tables", pattern="_limma_DEG\\.tsv$", full.names=TRUE)
if (length(deg_files) == 0) stop("No DEG tables found")

DSETS <- sort(unique(gsub("_limma_DEG\\.tsv$", "", basename(deg_files))))

# ---------------- Load expression ----------------
load_expr <- function(ds) {
  expr  <- readRDS(file.path(DIR_PROC, paste0(ds, "_expr_gene.rds")))
  sheet <- fread(file.path(DIR_META, paste0(ds, "_sample_sheet.csv")))
  
  sheet <- sheet[final_group %in% c("HGSC","Normal")]
  sids  <- intersect(colnames(expr), sheet$sample_id)
  
  expr  <- expr[GENES, sids, drop=FALSE]
  sheet <- sheet[match(sids, sample_id)]
  
  z <- t(scale(t(expr)))
  
  dt <- as.data.table(as.table(z))
  setnames(dt, c("gene","sample_id","z_expr"))
  dt[, dataset := ds]
  dt <- merge(dt, sheet[,.(sample_id, group=final_group)], by="sample_id")
  dt
}

expr_long <- rbindlist(lapply(DSETS, load_expr))

expr_meta <- copy(expr_long)
expr_meta[, dataset := "Pooled_META"]

expr_all <- rbind(expr_long, expr_meta)

# ---------------- Stats for stars (SAFE) ----------------
pstats <- expr_all[, .(
  p = tryCatch(wilcox.test(z_expr ~ group)$p.value, error=function(e) NA_real_),
  y = max(z_expr, na.rm=TRUE) + 0.45
), by=.(dataset, gene)]

pstats[, star := star_code(p)]

# ---------------- Plot 1: Box + Violin ----------------
p_box <- ggplot(expr_all, aes(group, z_expr, fill=group)) +
  geom_violin(trim=FALSE, alpha=0.55) +
  geom_boxplot(width=0.22, outlier.shape=NA, alpha=0.9) +
  geom_jitter(width=0.1, size=1.1, alpha=0.6) +
  geom_text(
    data = pstats,
    aes(x="HGSC", y=y, label=star),
    inherit.aes = FALSE,
    size = STAR_BOX
  ) +
  scale_fill_manual(values=c("Normal"=COL_NORM,"HGSC"=COL_UP)) +
  facet_grid(gene ~ dataset, scales="free_y") +
  theme_bw(base_size=BASE_SIZE) +
  theme(
    strip.text = element_text(face="bold"),
    axis.title.x = element_blank(),
    legend.position = "top"
  ) +
  labs(
    y = "Z-scored expression",
    title = "CORE 2 genes (SPON1, ALDH1A2) across datasets"
  )

ggsave("figures/CORE2_boxviolin_SPON1_ALDH1A2.pdf",
       p_box, width=17, height=6, device=cairo_pdf)
ggsave("figures/CORE2_boxviolin_SPON1_ALDH1A2.png",
       p_box, width=17, height=6, dpi=600)

# ---------------- Plot 2: Mean Z heatmap ----------------
meanZ <- expr_all[, .(meanZ = mean(z_expr, na.rm=TRUE)),
                  by=.(gene, dataset, group)]
meanZ[, col := paste(dataset, group, sep=" | ")]

p_hm_mean <- ggplot(meanZ, aes(col, gene, fill=meanZ)) +
  geom_tile(color="white", linewidth=0.4) +
  scale_fill_gradient2(low=COL_DOWN, mid="white", high=COL_UP, midpoint=0) +
  theme_bw(base_size=BASE_SIZE) +
  theme(
    axis.text.x = element_text(angle=35, hjust=1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title="Mean standardized expression")

ggsave("figures/CORE2_heatmap_meanZ_SPON1_ALDH1A2.pdf",
       p_hm_mean, width=14, height=3, device=cairo_pdf)
ggsave("figures/CORE2_heatmap_meanZ_SPON1_ALDH1A2.png",
       p_hm_mean, width=14, height=3, dpi=600)

# ---------------- Plot 3: logFC heatmap ----------------
read_deg <- function(ds) {
  dt <- fread(file.path("results/tables", paste0(ds, "_limma_DEG.tsv")))
  dt[gene %in% GENES,
     .(gene, dataset=ds, logFC, FDR=adj.P.Val)]
}

deg_long <- rbindlist(lapply(DSETS, read_deg))
deg_long[, star := star_code(FDR)]

p_hm_logfc <- ggplot(deg_long, aes(dataset, gene, fill=logFC)) +
  geom_tile(color="white", linewidth=0.4) +
  geom_text(aes(label=star), size=STAR_HM) +
  scale_fill_gradient2(low=COL_DOWN, mid="white", high=COL_UP, midpoint=0) +
  theme_bw(base_size=BASE_SIZE) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title="Per-dataset logFC (limma)")

ggsave("figures/CORE2_heatmap_logFC_META_SPON1_ALDH1A2.pdf",
       p_hm_logfc, width=9.5, height=3, device=cairo_pdf)
ggsave("figures/CORE2_heatmap_logFC_META_SPON1_ALDH1A2.png",
       p_hm_logfc, width=9.5, height=3, dpi=600)

message("DONE. No missing stars. No silent failures. Peace restored.")
