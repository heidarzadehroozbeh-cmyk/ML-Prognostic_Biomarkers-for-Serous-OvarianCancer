# scripts/05_meta_deg_consensus.R
# Stage 05: Meta-analysis across mRNA GEO datasets (GPL570) to derive consensus DEGs
#
# Inputs:
#   results/tables/<GSE>_limma_all.tsv  (from Stage 04)
#
# Outputs:
#   results/tables/META_all_genes.tsv
#   results/tables/META_consensus_STRICT.tsv
#   results/tables/META_consensus_RELAXED.tsv
#
# Notes:
#   - We meta-analyze study-level logFC using random-effects (REML) via metafor::rma.uni
#   - SE is approximated as |logFC|/|t| (from limma topTable)
#   - STRICT consensus: k==3 + same direction in all 3 + meta_FDR < DEG_FDR + |meta_logFC| > DEG_abs_logFC
#   - RELAXED consensus: k>=2 + majority direction + meta_FDR < DEG_FDR  (no strict effect-size cutoff)
#
# Run:
#   source("scripts/05_meta_deg_consensus.R")

source("scripts/01_params.R")

suppressPackageStartupMessages({
  library(data.table)
  library(metafor)
})

# ---------- helper: read limma table ----------
read_limma <- function(gse_id) {
  f <- file.path(DIR_TAB, paste0(gse_id, "_limma_all.tsv"))
  if (!file.exists(f)) stop("Missing lim hooking file: ", f)
  dt <- fread(f)
  df <- as.data.frame(dt)
  
  # standardize required columns
  req <- c("gene", "logFC", "t", "P.Value", "adj.P.Val")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop("Missing columns in ", f, ": ", paste(miss, collapse=", "))
  
  df
}

# ---------- helper: safe SE from limma output ----------
calc_se_from_t <- function(logFC, tstat) {
  # se ≈ |logFC|/|t|
  se <- abs(logFC) / abs(tstat)
  se[!is.finite(se)] <- NA
  se[se == 0] <- NA
  se
}

# ---------- load all studies ----------
studies <- GSE_mrna
tt_list <- lapply(studies, read_limma)
names(tt_list) <- studies

# union genes
all_genes <- Reduce(union, lapply(tt_list, \(x) x$gene))

message("Total unique genes across studies: ", length(all_genes))

# ---------- meta-analysis per gene ----------
meta_rows <- vector("list", length(all_genes))

for (i in seq_along(all_genes)) {
  g <- all_genes[i]
  
  per_study <- lapply(studies, function(sid) {
    tt <- tt_list[[sid]]
    r <- tt[tt$gene == g, ]
    if (nrow(r) == 0) return(NULL)
    
    se <- calc_se_from_t(r$logFC, r$t)
    data.frame(
      study = sid,
      logFC = as.numeric(r$logFC),
      se    = as.numeric(se),
      t     = as.numeric(r$t),
      p     = as.numeric(r$P.Value),
      fdr   = as.numeric(r$adj.P.Val),
      stringsAsFactors = FALSE
    )
  })
  
  per_study <- do.call(rbind, per_study)
  if (is.null(per_study) || nrow(per_study) < 2) {
    meta_rows[[i]] <- NULL
    next
  }
  
  # remove entries with missing SE
  per_study <- per_study[!is.na(per_study$se), , drop = FALSE]
  if (nrow(per_study) < 2) {
    meta_rows[[i]] <- NULL
    next
  }
  
  # random-effects meta
  m <- tryCatch(
    rma.uni(yi = per_study$logFC, sei = per_study$se, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(m)) {
    meta_rows[[i]] <- NULL
    next
  }
  
  # direction concordance
  dirs <- sign(per_study$logFC)
  dir_all_same <- (length(unique(dirs)) == 1)
  
  # majority direction (for k>=2)
  dir_majority <- ifelse(sum(dirs > 0) > sum(dirs < 0), 1,
                         ifelse(sum(dirs < 0) > sum(dirs > 0), -1, 0))
  
  # add per-study columns wide (optional but useful)
  wide <- setNames(as.list(rep(NA_real_, length(studies)*2)), c(
    as.vector(rbind(paste0("logFC_", studies), paste0("fdr_", studies)))
  ))
  for (sid in studies) {
    rr <- per_study[per_study$study == sid, ]
    if (nrow(rr) == 1) {
      wide[[paste0("logFC_", sid)]] <- rr$logFC
      wide[[paste0("fdr_", sid)]]   <- rr$fdr
    }
  }
  
  meta_rows[[i]] <- cbind(
    data.frame(
      gene = g,
      k = nrow(per_study),
      meta_logFC = as.numeric(m$b),
      meta_se    = as.numeric(m$se),
      meta_p     = as.numeric(m$pval),
      tau2       = as.numeric(m$tau2),
      I2         = as.numeric(m$I2),
      dir_all_same = dir_all_same,
      dir_majority = dir_majority,
      stringsAsFactors = FALSE
    ),
    as.data.frame(wide)
  )
}

meta <- do.call(rbind, meta_rows)
if (is.null(meta) || nrow(meta) == 0) stop("Meta table is empty; check Stage 04 outputs.")

meta$meta_FDR <- p.adjust(meta$meta_p, method = "BH")

# Save full meta table
out_all <- file.path(DIR_TAB, "META_all_genes.tsv")
write.table(meta, out_all, sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: ", out_all)

# ---------- CONSENSUS definitions ----------
# STRICT: all 3 studies available + same direction in all 3 + meta_FDR < cutoff + |meta_logFC| > cutoff
cons_strict <- subset(meta,
                      k == length(studies) &
                        dir_all_same == TRUE &
                        meta_FDR < DEG_FDR &
                        abs(meta_logFC) > DEG_abs_logFC)

out_strict <- file.path(DIR_TAB, "META_consensus_STRICT.tsv")
write.table(cons_strict, out_strict, sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: ", out_strict, " | N_STRICT = ", nrow(cons_strict))

# RELAXED: k>=2 + majority direction exists (dir_majority != 0) + meta_FDR < cutoff
# (informative)
cons_relaxed <- subset(meta,
                       k >= 2 &
                         dir_majority != 0 &
                         meta_FDR < DEG_FDR)

out_relaxed <- file.path(DIR_TAB, "META_consensus_RELAXED.tsv")
write.table(cons_relaxed, out_relaxed, sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: ", out_relaxed, " | N_RELAXED = ", nrow(cons_relaxed))

# optional: quick up/down lists from STRICT
if (nrow(cons_strict) > 0) {
  up <- cons_strict$gene[cons_strict$meta_logFC > 0]
  dn <- cons_strict$gene[cons_strict$meta_logFC < 0]
  write.table(up, file.path(DIR_TAB, "META_STRICT_up_genes.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(dn, file.path(DIR_TAB, "META_STRICT_down_genes.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  message("Wrote STRICT up/down gene lists.")
}

message("\nStage 05 complete.")
