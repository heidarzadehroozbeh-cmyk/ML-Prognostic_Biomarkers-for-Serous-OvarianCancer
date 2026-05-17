# scripts/21_volcano_colored_limma_Q1.R
# Colored volcano plots for limma tables (CRAN-only)
# Fix: ensure NS (grey) points appear by:
#  1) auto-detecting if input is filtered; try to switch to unfiltered limma table if available
#  2) plotting in layers: NS first, then Down, then Up (so NS not hidden)
#  3) handling zeros in p/FDR (avoid dropping)
# Works with tables containing logFC + adj.P.Val/FDR (preferred) or p-value (fallback)


source("scripts/01_params.R")

options(timeout = 600)
pkgs <- c("data.table","ggplot2","stringr","patchwork")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
if (!requireNamespace("ggrepel", quietly=TRUE)) install.packages("ggrepel")

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(stringr); library(patchwork); library(ggrepel)
})

# ---------------- USER SETTINGS ----------------
FDR_CUT <- 0.05
ABS_LOGFC_CUT <- 1.0
N_LABELS <- 12
LABEL_ONLY_SIGNIFICANT <- TRUE

# If NS proportion is too small, try switching to an unfiltered limma file (if exists)
NS_MIN_PROP <- 0.02
TRY_FIND_UNFILTERED <- TRUE

COL_UP   <- "#E64B35"   # warm
COL_DOWN <- "#00BFC4"   # turquoise
COL_NS   <- "#BDBDBD"   # grey

# ---------------- Helpers ----------------
pick_col <- function(dt, candidates) {
  cand <- candidates[candidates %in% names(dt)]
  if (length(cand) == 0) return(NA_character_)
  cand[1]
}

safe_tag <- function(x) {
  x <- gsub("<=|≥|≤", "le", x)
  x <- gsub("<", "lt", x)
  x <- gsub(">", "gt", x)
  x <- gsub("=", "eq", x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x
}

pdf_dev <- function() {
  if (capabilities("cairo")) return(cairo_pdf)
  return(pdf)
}

# Find an unfiltered limma table for same dataset prefix
find_unfiltered <- function(infile) {
  base <- basename(infile)
  ds <- gsub("_limma_DEG\\.tsv$", "", base)
  
  # candidate patterns (same directory)
  dir0 <- dirname(infile)
  cands <- c(
    file.path(dir0, paste0(ds, "_limma_ALL.tsv")),
    file.path(dir0, paste0(ds, "_limma_all.tsv")),
    file.path(dir0, paste0(ds, "_limma_full.tsv")),
    file.path(dir0, paste0(ds, "_limma_FULL.tsv")),
    file.path(dir0, paste0(ds, "_limma_table.tsv")),
    file.path(dir0, paste0(ds, "_limma_results.tsv")),
    file.path(dir0, paste0(ds, "_limma.tsv"))
  )
  cands <- cands[file.exists(cands)]
  if (length(cands) > 0) return(cands[1])
  
  # broader search: any limma*.tsv with same ds but NOT containing "_DEG"
  all_hits <- list.files(dir0, pattern = paste0("^", ds, "_limma.*\\.tsv$"), full.names = TRUE)
  all_hits <- all_hits[!grepl("_DEG\\.tsv$", all_hits)]
  if (length(all_hits) > 0) return(all_hits[1])
  
  return(NA_character_)
}

# Convert p-like values to numeric safely and handle zeros
sanitize_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  TRUE
}

fix_zero_p <- function(p) {
  p <- as.numeric(p)
  p[!is.finite(p)] <- NA_real_
  # replace <=0 with a small positive (based on smallest positive)
  min_pos <- suppressWarnings(min(p[p > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-300
  p[p <= 0] <- min_pos * 1e-1
  p
}

# ---------------- Locate input DEG files ----------------
deg_files <- list.files("results/tables", pattern="_limma_DEG\\.tsv$", full.names=TRUE)
if (length(deg_files) == 0) stop("No DEG files found in results/tables/*_limma_DEG.tsv")
message("[INFO] Found DEG files:\n", paste(deg_files, collapse="\n"))

# Optional: FINAL25 list for highlighting
final25_files <- list.files("figures", pattern="FINAL25_gene_list_.*\\.tsv$", full.names=TRUE)
final25 <- NULL
if (length(final25_files) > 0) {
  message("[INFO] Found FINAL25 list: ", final25_files[1])
  final25_dt <- fread(final25_files[1])
  gcol_f <- pick_col(final25_dt, c("gene","Gene","symbol","SYMBOL","GeneSymbol","hgnc_symbol"))
  if (is.na(gcol_f)) gcol_f <- names(final25_dt)[1]
  final25 <- unique(as.character(final25_dt[[gcol_f]]))
} else {
  message("[INFO] FINAL25 list not found in figures/. Proceeding without highlighting.")
}

# ---------------- Output dir ----------------
DIR_FIG <- if (exists("DIR_FIG")) DIR_FIG else "figures"
out_dir <- file.path(DIR_FIG, "VOLCANO_Q1")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# ---------------- Volcano function ----------------
make_volcano <- function(infile) {
  
  load_one <- function(f) {
    dt <- fread(f)
    
    gene_col  <- pick_col(dt, c("gene","Gene","symbol","SYMBOL","GeneSymbol","hgnc_symbol"))
    logfc_col <- pick_col(dt, c("meta_logFC","logFC","log2FC","log2FoldChange","effect","Estimate","beta"))
    fdr_col   <- pick_col(dt, c("meta_fdr","FDR","padj","adj.P.Val","adj.Pval","qvalue","q_value","BH_FDR"))
    p_col     <- pick_col(dt, c("P.Value","pvalue","p_value","pval","P"))
    
    if (is.na(gene_col)) gene_col <- names(dt)[1]
    if (is.na(logfc_col)) stop("Cannot find logFC column in: ", f)
    
    # choose p for plotting: prefer adjusted if available
    p_used_col <- if (!is.na(fdr_col)) fdr_col else p_col
    if (is.na(p_used_col)) stop("Cannot find FDR/adj.P.Val nor p-value column in: ", f)
    
    dt[, gene := as.character(get(gene_col))]
    dt[, logFC := suppressWarnings(as.numeric(get(logfc_col)))]
    
    dt[, p_used := suppressWarnings(as.numeric(get(p_used_col)))]
    dt <- dt[!is.na(gene) & gene != "" & is.finite(logFC) & is.finite(p_used)]
    dt[, p_used := fix_zero_p(p_used)]
    dt <- dt[p_used > 0]
    
    # volcano y
    dt[, neglog10 := -log10(p_used)]
    
    # status
    dt[, status := "NS"]
    dt[p_used <= FDR_CUT & logFC >=  ABS_LOGFC_CUT, status := "Up"]
    dt[p_used <= FDR_CUT & logFC <= -ABS_LOGFC_CUT, status := "Down"]
    dt[, status := factor(status, levels=c("NS","Down","Up"))]
    
    # final25 highlight
    dt[, is_final := 0L]
    if (!is.null(final25)) dt[gene %in% final25, is_final := 1L]
    
    list(dt=dt, pcol=p_used_col, used_adj = !is.na(fdr_col))
  }
  
  base <- basename(infile)
  ds <- gsub("_limma_DEG\\.tsv$", "", base)
  
  obj <- load_one(infile)
  dt <- obj$dt
  p_used_col <- obj$pcol
  used_adj <- obj$used_adj
  
  # If NS is almost absent, try to switch to unfiltered limma table
  prop_ns <- mean(dt$status == "NS", na.rm = TRUE)
  if (TRY_FIND_UNFILTERED && is.finite(prop_ns) && prop_ns < NS_MIN_PROP) {
    alt <- find_unfiltered(infile)
    if (!is.na(alt)) {
      message("[INFO] ", ds, ": NS proportion is low (", sprintf("%.3f", prop_ns),
              "). Switching to unfiltered file: ", alt)
      obj2 <- load_one(alt)
      dt2 <- obj2$dt
      prop_ns2 <- mean(dt2$status == "NS", na.rm = TRUE)
      # use alt if it improves NS presence or if dt2 is much larger
      if (nrow(dt2) > nrow(dt) || prop_ns2 > prop_ns) {
        dt <- dt2
        p_used_col <- obj2$pcol
        used_adj <- obj2$used_adj
        prop_ns <- prop_ns2
      }
    } else {
      message("[WARN] ", ds, ": NS proportion is low and no unfiltered limma table found. ",
              "This usually means your *_limma_DEG.tsv is already filtered to significant genes only.")
    }
  }
  
  # label selection
  lab_dt <- copy(dt)
  if (LABEL_ONLY_SIGNIFICANT) lab_dt <- lab_dt[status != "NS"]
  lab_dt[, abs_logFC := abs(logFC)]
  
  if (!is.null(final25) && any(lab_dt$is_final==1L)) {
    lab_dt <- lab_dt[is_final==1L]
    setorder(lab_dt, -neglog10, -abs_logFC)
    lab_dt <- unique(lab_dt, by="gene")[1:min(N_LABELS, .N)]
  } else {
    setorder(lab_dt, -neglog10, -abs_logFC)
    lab_dt <- unique(lab_dt, by="gene")[1:min(N_LABELS, .N)]
  }
  
  # axis label depends on p used
  y_lab <- if (used_adj) {
    expression(-log[10]("adj.P.Val / FDR"))
  } else {
    expression(-log[10]("p-value"))
  }
  
  # IMPORTANT: draw in layers so NS is visible
  p <- ggplot() +
    geom_point(
      data = dt[status=="NS"],
      aes(x=logFC, y=neglog10),
      color = COL_NS, alpha=0.60, size=1.10
    ) +
    geom_point(
      data = dt[status=="Down"],
      aes(x=logFC, y=neglog10),
      color = COL_DOWN, alpha=0.80, size=1.25
    ) +
    geom_point(
      data = dt[status=="Up"],
      aes(x=logFC, y=neglog10),
      color = COL_UP, alpha=0.80, size=1.25
    ) +
    { if (!is.null(final25)) geom_point(
      data = dt[is_final==1L],
      aes(x=logFC, y=neglog10),
      shape=21, fill=NA, color="black", stroke=0.55, size=1.90, alpha=0.95
    ) } +
    geom_vline(xintercept=c(-ABS_LOGFC_CUT, ABS_LOGFC_CUT),
               linetype="dashed", linewidth=0.5, color="grey40") +
    geom_hline(yintercept=-log10(FDR_CUT),
               linetype="dashed", linewidth=0.5, color="grey40") +
    ggrepel::geom_text_repel(
      data=lab_dt,
      aes(x=logFC, y=neglog10, label=gene),
      size=3.1,
      max.overlaps=Inf,
      box.padding=0.35,
      point.padding=0.25,
      segment.color="grey55",
      segment.size=0.3
    ) +
    theme_bw(base_size=12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size=14, face="bold"),
      plot.subtitle = element_text(size=11),
      legend.position = "none",
      plot.margin = margin(10, 12, 10, 10)
    ) +
    labs(
      title = paste0(ds, " | Volcano (limma)"),
      subtitle = sprintf(
        "%s ≤ %.2g and |logFC| ≥ %.2f | Up=warm, Down=turquoise, NS=grey | NS%%=%.1f",
        if (used_adj) "BH-FDR" else "p-value",
        FDR_CUT, ABS_LOGFC_CUT, 100*prop_ns
      ),
      x = "log2 fold change",
      y = y_lab
    )
  
  tag2 <- safe_tag(sprintf("CUT_%.2g_%.2f", FDR_CUT, ABS_LOGFC_CUT))
  out_pdf <- file.path(out_dir, paste0(ds, "_VOLCANO_", tag2, ".pdf"))
  out_png <- file.path(out_dir, paste0(ds, "_VOLCANO_", tag2, ".png"))
  
  ggsave(out_pdf, p, width=7.6, height=6.4, units="in", device=pdf_dev(), limitsize=FALSE)
  ggsave(out_png, p, width=7.6, height=6.4, units="in", dpi=600, bg="white", limitsize=FALSE)
  
  message("[OK] ", ds, " | NS proportion: ", sprintf("%.3f", prop_ns))
  message("[OK] Saved: ", out_pdf)
  message("[OK] Saved: ", out_png)
  
  list(ds=ds, plot=p)
}

plots <- lapply(deg_files, make_volcano)
p_list <- lapply(plots, `[[`, "plot")
names(p_list) <- vapply(plots, `[[`, character(1), "ds")

p_comb <- wrap_plots(p_list, ncol=length(p_list)) +
  plot_annotation(
    title = "Volcano plots across cohorts (limma)",
    subtitle = "Colored by directionality and significance; grey = NS. FINAL genes (if provided) are outlined in black."
  )

tag_all <- safe_tag(sprintf("CUT_%.2g_%.2f", FDR_CUT, ABS_LOGFC_CUT))
out_pdf_all <- file.path(out_dir, paste0("ALL_VOLCANO_PANEL_", tag_all, ".pdf"))
out_png_all <- file.path(out_dir, paste0("ALL_VOLCANO_PANEL_", tag_all, ".png"))

ggsave(out_pdf_all, p_comb, width=20.0, height=6.6, units="in", device=pdf_dev(), limitsize=FALSE)
ggsave(out_png_all, p_comb, width=20.0, height=6.6, units="in", dpi=600, bg="white", limitsize=FALSE)

message("DONE. Combined panel saved to: ", out_pdf_all)
