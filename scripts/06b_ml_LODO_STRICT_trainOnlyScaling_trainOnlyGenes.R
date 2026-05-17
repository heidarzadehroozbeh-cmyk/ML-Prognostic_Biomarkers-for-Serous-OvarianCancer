# scripts/06b_ml_LODO_STRICT_trainOnlyScaling_trainOnlyGenes.R
source("scripts/01_params.R")

pkgs <- c("data.table","glmnet","pROC","Matrix","matrixStats")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(pROC); library(Matrix); library(matrixStats)
})

# ---------------- helpers ----------------
read_limma_table <- function(ds) {
  f <- file.path(DIR_TAB, paste0(ds, "_limma_all.tsv"))
  if (!file.exists(f)) stop("Missing limma table: ", f)
  tt <- fread(f)
  
  gene_col  <- if ("gene"  %in% names(tt)) "gene"  else stop("No gene column in ", f)
  logfc_col <- if ("logFC" %in% names(tt)) "logFC" else stop("No logFC column in ", f)
  
  fdr_candidates <- c("FDR","adj.P.Val","adjP","adj_p","fdr","padj")
  fdr_col <- fdr_candidates[fdr_candidates %in% names(tt)]
  if (length(fdr_col)==0) stop("No FDR/adj.P.Val column found in ", f)
  fdr_col <- fdr_col[1]
  
  out <- tt[, .(
    gene = as.character(get(gene_col)),
    logFC = as.numeric(get(logfc_col)),
    fdr = as.numeric(get(fdr_col))
  )]
  out <- out[is.finite(logFC) & is.finite(fdr)]
  out
}

# TRAIN-only candidate genes (no leakage), with safe unique column names
candidate_genes_train_only <- function(train_ds, fdr_cut=0.05, min_abs_logFC=0, mode=c("strict","relaxed")) {
  mode <- match.arg(mode)
  tabs <- lapply(train_ds, function(ds) {
    tt <- read_limma_table(ds)
    setnames(tt, c("logFC","fdr"), c(paste0("logFC_", ds), paste0("fdr_", ds)))
    tt
  })
  
  common <- Reduce(intersect, lapply(tabs, `[[`, "gene"))
  tabs <- lapply(tabs, function(x) x[gene %in% common])
  
  m <- Reduce(function(a,b) merge(a,b, by="gene", all=FALSE), tabs)
  
  logfc_cols <- paste0("logFC_", train_ds)
  fdr_cols   <- paste0("fdr_",   train_ds)
  
  # filter per dataset
  keep <- rep(TRUE, nrow(m))
  for (j in seq_along(train_ds)) {
    keep <- keep & (m[[fdr_cols[j]]] < fdr_cut) & (abs(m[[logfc_cols[j]]]) >= min_abs_logFC)
  }
  if (!any(keep)) return(character(0))
  mm <- m[keep, ]
  
  # direction concordance
  dirs <- sign(as.matrix(mm[, ..logfc_cols]))
  same_dir <- apply(dirs, 1, function(v) length(unique(v))==1)
  mm <- mm[same_dir, ]
  
  mm$gene
}

load_expr_and_labels <- function(ds) {
  expr_path  <- file.path(DIR_PROC, paste0(ds, "_expr_gene.rds"))
  sheet_path <- file.path(DIR_META, paste0(ds, "_sample_sheet.csv"))
  if (!file.exists(expr_path))  stop("Missing expr: ", expr_path)
  if (!file.exists(sheet_path)) stop("Missing sheet: ", sheet_path)
  
  expr <- readRDS(expr_path)  # genes x samples
  sheet <- read.csv(sheet_path)
  
  if (!("final_group" %in% names(sheet))) stop("final_group missing in ", sheet_path)
  if (!("sample_id" %in% names(sheet))) stop("sample_id missing in ", sheet_path)
  
  sheet$final_group[sheet$final_group==""] <- NA
  sheet <- sheet[!is.na(sheet$final_group) & sheet$final_group %in% c("HGSC","Normal"), ]
  
  sids <- intersect(colnames(expr), sheet$sample_id)
  if (length(sids) < 6) stop("Too few matched samples in ", ds)
  
  sheet <- sheet[match(sids, sheet$sample_id), ]
  expr  <- expr[, sids, drop=FALSE]
  
  y <- ifelse(sheet$final_group=="HGSC", 1L, 0L)
  list(expr=expr, y=y, sample_id=sids)
}

scale_train_apply_test <- function(Xtr, Xte) {
  mu <- Matrix::colMeans(Xtr)
  sd <- sqrt(Matrix::colMeans((Xtr - Matrix(mu, nrow=nrow(Xtr), ncol=ncol(Xtr), byrow=TRUE))^2))
  sd[sd==0 | !is.finite(sd)] <- 1
  
  Xtr_s <- (Xtr - Matrix(mu, nrow=nrow(Xtr), ncol=ncol(Xtr), byrow=TRUE)) /
    Matrix(sd, nrow=nrow(Xtr), ncol=ncol(Xtr), byrow=TRUE)
  
  Xte_s <- (Xte - Matrix(mu, nrow=nrow(Xte), ncol=ncol(Xte), byrow=TRUE)) /
    Matrix(sd, nrow=nrow(Xte), ncol=ncol(Xte), byrow=TRUE)
  
  list(Xtr=Xtr_s, Xte=Xte_s)
}

# ---------------- settings ----------------
alpha <- 0.5
fdr_cut_strict <- 0.05
fdr_cut_relaxed <- 0.10
min_abs_logFC <- 0

stability_B <- 80
subsample <- 0.80
freq_cut_main <- 0.60

MIN_PANEL <- 10           # minimum genes to keep as a usable panel
TOPK_PANEL <- 25          # if frequency threshold yields small panel, take TopK by frequency
TOP_VAR_GENES <- 3000     # unsupervised train-only variance filter (helps stability)

set.seed(1)

# ---------------- datasets ----------------
datasets <- GSE_mrna
datasets <- datasets[
  file.exists(file.path(DIR_PROC, paste0(datasets, "_expr_gene.rds"))) &
    file.exists(file.path(DIR_META, paste0(datasets, "_sample_sheet.csv"))) &
    file.exists(file.path(DIR_TAB,  paste0(datasets, "_limma_all.tsv")))
]
if (length(datasets) < 3) stop("Need >=3 datasets with expr+sheet+limma_all.tsv.")

DL <- lapply(datasets, load_expr_and_labels)
names(DL) <- datasets

# ---------------- LODO loop ----------------
pred_rows  <- list()
auc_rows   <- list()
panel_rows <- list()

for (test_ds in datasets) {
  train_ds <- setdiff(datasets, test_ds)
  
  # (A) TRAIN-only candidate genes (strict then fallback relaxed) - still no leakage
  cand <- candidate_genes_train_only(train_ds, fdr_cut=fdr_cut_strict, min_abs_logFC=min_abs_logFC, mode="strict")
  if (length(cand) < 200) {
    message("[INFO] Few strict candidates in test=", test_ds, " (n=", length(cand), "), trying relaxed fdr<0.10")
    cand2 <- candidate_genes_train_only(train_ds, fdr_cut=fdr_cut_relaxed, min_abs_logFC=min_abs_logFC, mode="relaxed")
    if (length(cand2) > length(cand)) cand <- cand2
  }
  
  # (B) genes available across train+test
  common_all <- Reduce(intersect, lapply(c(train_ds, test_ds), function(ds) rownames(DL[[ds]]$expr)))
  
  genes_use <- common_all
  if (length(cand) >= 200) genes_use <- intersect(common_all, cand)
  
  if (length(genes_use) < 200) {
    message("[WARN] genes_use small (", length(genes_use), ") in test=", test_ds, " ; using common_all instead.")
    genes_use <- common_all
  }
  
  # (C) build train/test matrices (samples x genes)
  Xtr0 <- do.call(rbind, lapply(train_ds, function(ds) t(DL[[ds]]$expr[genes_use, , drop=FALSE])))
  ytr  <- unlist(lapply(train_ds, function(ds) DL[[ds]]$y))
  
  Xte0 <- t(DL[[test_ds]]$expr[genes_use, , drop=FALSE])
  yte  <- DL[[test_ds]]$y
  
  # (D) unsupervised train-only variance filter (no outcome leakage)
  if (ncol(Xtr0) > TOP_VAR_GENES) {
    v <- matrixStats::colVars(as.matrix(Xtr0))
    keep_idx <- order(v, decreasing=TRUE)[1:TOP_VAR_GENES]
    Xtr0 <- Xtr0[, keep_idx, drop=FALSE]
    Xte0 <- Xte0[, keep_idx, drop=FALSE]
    genes_use <- colnames(Xtr0)
    message("[INFO] Variance filter applied: kept ", length(genes_use), " genes in test=", test_ds)
  }
  
  Xtr <- Matrix(Xtr0, sparse=FALSE)
  Xte <- Matrix(Xte0, sparse=FALSE)
  
  # (E) TRAIN-only scaling (no test leakage)
  sc <- scale_train_apply_test(Xtr, Xte)
  Xtr <- sc$Xtr; Xte <- sc$Xte
  
  # (F) CV on TRAIN only: try lambda.1se then fallback to lambda.min if too sparse
  cvfit <- cv.glmnet(Xtr, ytr, family="binomial", alpha=alpha, nfolds=5, type.measure="auc")
  lam_1se <- cvfit$lambda.1se
  lam_min <- cvfit$lambda.min
  
  choose_lambda <- function(lam) {
    # stability selection on TRAIN only
    sel_count <- setNames(rep(0L, ncol(Xtr)), colnames(Xtr))
    ntr <- nrow(Xtr)
    for (b in 1:stability_B) {
      idx <- sample.int(ntr, size=floor(subsample*ntr), replace=FALSE)
      fit_b <- glmnet(Xtr[idx, , drop=FALSE], ytr[idx], family="binomial", alpha=alpha, lambda=lam)
      beta <- as.matrix(coef(fit_b))[-1,1]
      picked <- names(beta)[beta != 0]
      sel_count[picked] <- sel_count[picked] + 1L
    }
    freq <- sel_count / stability_B
    list(freq=freq)
  }
  
  ss1 <- choose_lambda(lam_1se)
  freq <- ss1$freq
  
  panel <- names(freq)[freq >= freq_cut_main]
  if (length(panel) < MIN_PANEL) {
    # fallback: take top genes by freq
    ord <- order(freq, decreasing=TRUE)
    panel <- names(freq)[ord][1:min(TOPK_PANEL, length(ord))]
    panel <- panel[is.finite(freq[panel])]
    message("[WARN] threshold panel<", MIN_PANEL, " in test=", test_ds, " ; using Top-", length(panel), " by stability freq (lambda.1se)")
  }
  
  # if still too small/non-informative, re-run with lambda.min
  if (length(panel) < MIN_PANEL) {
    message("[INFO] Trying lambda.min in test=", test_ds)
    ss2 <- choose_lambda(lam_min)
    freq <- ss2$freq
    ord <- order(freq, decreasing=TRUE)
    panel <- names(freq)[ord][1:min(TOPK_PANEL, length(ord))]
    panel <- panel[is.finite(freq[panel])]
    message("[WARN] Using Top-", length(panel), " by stability freq (lambda.min) in test=", test_ds)
  }
  
  if (length(panel) < 5) stop("Panel still too small in test=", test_ds, " (", length(panel), ")")
  
  # (G) refit on TRAIN with panel only, predict TEST
  Xtr_p <- Xtr[, panel, drop=FALSE]
  Xte_p <- Xte[, panel, drop=FALSE]
  
  # choose lambda used for final fit: use lambda.1se (more regularized); acceptable because panel already constrained
  final <- glmnet(Xtr_p, ytr, family="binomial", alpha=alpha, lambda=lam_1se)
  prob  <- as.numeric(predict(final, Xte_p, type="response"))
  auc   <- as.numeric(pROC::auc(pROC::roc(yte, prob, quiet=TRUE)))
  
  pred_rows[[test_ds]] <- data.frame(
    test_dataset=test_ds,
    sample_id=DL[[test_ds]]$sample_id,
    y=yte,
    prob=prob,
    stringsAsFactors=FALSE
  )
  auc_rows[[test_ds]] <- data.frame(test_dataset=test_ds, AUC=auc, panel_size=length(panel))
  panel_rows[[test_ds]] <- data.frame(test_dataset=test_ds, gene=panel, freq=freq[panel], stringsAsFactors=FALSE)
  
  message("STRICT LODO | test=", test_ds, " | AUC=", signif(auc,4), " | panel=", length(panel))
}

pred_out  <- rbindlist(pred_rows)
auc_out   <- rbindlist(auc_rows)
panel_out <- rbindlist(panel_rows)

dir.create(DIR_MOD, showWarnings=FALSE, recursive=TRUE)
fwrite(pred_out,  file.path(DIR_MOD, "LODO_STRICT_predictions.tsv"), sep="\t")
fwrite(auc_out,   file.path(DIR_MOD, "LODO_STRICT_AUC.csv"))
fwrite(panel_out, file.path(DIR_MOD, "LODO_STRICT_panels.tsv"), sep="\t")

message("\nSaved:")
message(" - results/models/LODO_STRICT_predictions.tsv")
message(" - results/models/LODO_STRICT_AUC.csv")
message(" - results/models/LODO_STRICT_panels.tsv")
