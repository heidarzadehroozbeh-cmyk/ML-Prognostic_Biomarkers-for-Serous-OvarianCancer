# scripts/06_ml_diagnostic_LODO.R
# Stage 06: Diagnostic ML (HGSC vs Normal) with cross-study LODO validation
# - Baseline: Elastic Net on candidate gene set
# - Nested stability selection: select a smaller gene panel within each LODO train set

source("scripts/01_params.R")

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(data.table)
  library(matrixStats)
  library(doParallel)
  library(foreach)
})

# ---------------------------
# 0) Load candidate genes
# ---------------------------
f_relaxed <- file.path(DIR_TAB, "META_consensus_RELAXED.tsv")
f_strict  <- file.path(DIR_TAB, "META_consensus_STRICT.tsv")

if (file.exists(f_relaxed)) {
  cand <- fread(f_relaxed)
  candidate_genes <- unique(cand$gene)
  message("Using RELAXED candidate genes: N = ", length(candidate_genes))
} else if (file.exists(f_strict)) {
  cand <- fread(f_strict)
  candidate_genes <- unique(cand$gene)
  message("RELAXED not found. Using STRICT candidate genes: N = ", length(candidate_genes))
} else {
  stop("No consensus gene file found. Run Stage 05 first.")
}

# ---------------------------
# 1) Build merged dataset (3 GEO mRNA sets)
# ---------------------------
load_one_dataset <- function(gse_id) {
  expr_gene <- readRDS(file.path(DIR_PROC, paste0(gse_id, "_expr_gene.rds")))  # genes x samples
  sheet <- read.csv(file.path(DIR_META, paste0(gse_id, "_sample_sheet.csv")))
  sheet$final_group[sheet$final_group == ""] <- NA
  sheet <- sheet[!is.na(sheet$final_group) & sheet$final_group %in% c("HGSC","Normal"), ]
  
  # subset to candidate genes present
  common <- intersect(candidate_genes, rownames(expr_gene))
  expr_sub <- expr_gene[common, sheet$sample_id, drop = FALSE]
  
  list(
    X = t(expr_sub),                         # samples x genes
    y = sheet$final_group,                   # HGSC/Normal
    ds = rep(gse_id, nrow(t(expr_sub)))
  )
}

dat_list <- lapply(GSE_mrna, load_one_dataset)
names(dat_list) <- GSE_mrna

# intersect genes across datasets to ensure same feature space
gene_sets <- lapply(dat_list, function(z) colnames(z$X))
genes_use <- Reduce(intersect, gene_sets)

if (length(genes_use) < 20) {
  stop("Too few overlapping genes across datasets after filtering: ", length(genes_use),
       "\nConsider loosening candidate gene criteria (use RELAXED) or check preprocessing.")
}

message("Final feature space (intersection across datasets): ", length(genes_use), " genes")

X <- do.call(rbind, lapply(dat_list, \(z) z$X[, genes_use, drop=FALSE]))
y_chr <- unlist(lapply(dat_list, \(z) z$y))
ds <- unlist(lapply(dat_list, \(z) z$ds))

y <- ifelse(y_chr == "HGSC", 1, 0)

# ---------------------------
# 2) Cross-study normalization: z-score within each dataset (gene-wise)
# ---------------------------
zscore_within_dataset <- function(X, ds_vec) {
  Xz <- X
  for (d in unique(ds_vec)) {
    idx <- which(ds_vec == d)
    Xz[idx, ] <- scale(X[idx, , drop = FALSE])
  }
  Xz
}
Xz <- zscore_within_dataset(X, ds)

# ---------------------------
# 3) Baseline LODO with Elastic Net
# ---------------------------
fit_predict_glmnet <- function(xtr, ytr, xte, alpha = 0.5) {
  cv <- cv.glmnet(xtr, ytr, family="binomial", alpha=alpha, nfolds=5)
  fit <- glmnet(xtr, ytr, family="binomial", alpha=alpha, lambda=cv$lambda.1se)
  p <- as.numeric(predict(fit, newx=xte, type="response"))
  list(prob=p, fit=fit, lambda=cv$lambda.1se)
}

lodo_baseline <- function(X, y, ds, alpha=0.5) {
  out <- list()
  perf <- data.frame(test_dataset=character(), AUC=numeric(), stringsAsFactors=FALSE)
  
  for (test_ds in unique(ds)) {
    train_idx <- which(ds != test_ds)
    test_idx  <- which(ds == test_ds)
    
    res <- fit_predict_glmnet(X[train_idx, , drop=FALSE], y[train_idx],
                              X[test_idx,  , drop=FALSE], alpha=alpha)
    
    roc_obj <- pROC::roc(y[test_idx], res$prob, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    
    perf <- rbind(perf, data.frame(test_dataset=test_ds, AUC=auc_val))
    out[[test_ds]] <- list(prob=res$prob, y=y[test_idx], lambda=res$lambda)
    message("Baseline LODO | Test=", test_ds, " | AUC=", round(auc_val, 4))
  }
  
  list(perf=perf, details=out)
}

baseline <- lodo_baseline(Xz, y, ds, alpha=0.5)

dir.create(DIR_MOD, showWarnings = FALSE, recursive = TRUE)
write.csv(baseline$perf, file.path(DIR_MOD, "LODO_baseline_AUC.csv"), row.names=FALSE)

message("Baseline mean AUC = ", round(mean(baseline$perf$AUC), 4))

# ---------------------------
# 4) Nested Stability Selection within each LODO fold
#    (panel selection only on training sets -> no leakage)
# ---------------------------
stability_select <- function(xtr, ytr, B=100, alpha=0.5, subsample=0.8, cores=4) {
  # returns selection frequency per gene
  n <- nrow(xtr); p <- ncol(xtr)
  genes <- colnames(xtr)
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  freq <- foreach(b = 1:B, .combine = `+`, .packages="glmnet") %dopar% {
    set.seed(1000 + b)
    idx <- sample(seq_len(n), size = ceiling(subsample*n), replace = FALSE)
    cv <- cv.glmnet(xtr[idx, , drop=FALSE], ytr[idx], family="binomial", alpha=alpha, nfolds=5)
    fit <- glmnet(xtr[idx, , drop=FALSE], ytr[idx], family="binomial", alpha=alpha, lambda=cv$lambda.1se)
    coef_vec <- as.matrix(coef(fit))[-1, 1]  # drop intercept
    as.numeric(coef_vec != 0)
  }
  
  stopCluster(cl)
  
  freq <- freq / B
  names(freq) <- genes
  freq
}

nested_lodo_with_panel <- function(X, y, ds, alpha=0.5,
                                   B=100, subsample=0.8, freq_cut=0.6, cores=4) {
  
  perf <- data.frame(test_dataset=character(), AUC=numeric(), n_panel=integer(), stringsAsFactors=FALSE)
  
  panel_list <- list()
  
  for (test_ds in unique(ds)) {
    message("\n--- Nested LODO fold: test = ", test_ds, " ---")
    train_idx <- which(ds != test_ds)
    test_idx  <- which(ds == test_ds)
    
    xtr <- X[train_idx, , drop=FALSE]
    ytr <- y[train_idx]
    xte <- X[test_idx,  , drop=FALSE]
    yte <- y[test_idx]
    
    # stability selection on TRAIN only
    freq <- stability_select(xtr, ytr, B=B, alpha=alpha, subsample=subsample, cores=cores)
    panel <- names(freq)[which(freq >= freq_cut)]
    
    # fallback: if too small, take top-N
    if (length(panel) < 5) {
      ord <- order(freq, decreasing = TRUE)
      panel <- names(freq)[ord[1:min(20, length(freq))]]
      message("Panel too small at cut; fallback to top-", length(panel))
    }
    
    # refit model on TRAIN using selected panel only
    res <- fit_predict_glmnet(xtr[, panel, drop=FALSE], ytr, xte[, panel, drop=FALSE], alpha=alpha)
    roc_obj <- pROC::roc(yte, res$prob, quiet=TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    
    perf <- rbind(perf, data.frame(test_dataset=test_ds, AUC=auc_val, n_panel=length(panel)))
    panel_list[[test_ds]] <- data.frame(gene=panel, freq=freq[panel], test_dataset=test_ds)
    
    message("Nested LODO | Test=", test_ds,
            " | AUC=", round(auc_val,4),
            " | Panel size=", length(panel))
  }
  
  list(perf=perf, panels=panel_list)
}

# Practical settings (if the system is more powerful, you can increase cores/B)
nested <- nested_lodo_with_panel(
  Xz, y, ds,
  alpha = 0.5,
  B = 80,           # Start with 80, can increase to 200 later
  subsample = 0.8,
  freq_cut = 0.6,
  cores = 4
)


write.csv(nested$perf, file.path(DIR_MOD, "LODO_nested_stability_AUC.csv"), row.names=FALSE)

# save per-fold panels
panel_df <- do.call(rbind, nested$panels)
write.table(panel_df, file.path(DIR_MOD, "LODO_nested_panels.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

message("\nStage 06 complete.")
message("Baseline AUC file: results/models/LODO_baseline_AUC.csv")
message("Nested stability AUC file: results/models/LODO_nested_stability_AUC.csv")
message("Panels file: results/models/LODO_nested_panels.tsv")
