# scripts/08_multimodel_gene_ranking_STRICT_LODO_top500.R
# Multi-model benchmarking + gene ranking on TRAIN-only Top-500 significant genes within each LODO fold
# Models: ElasticNet (glmnet), Linear SVM (e1071), RandomForest (ranger), XGBoost (xgboost), GBM (gbm)

source("scripts/01_params.R")

pkgs <- c("data.table","pROC","glmnet","Matrix","matrixStats","e1071","ranger","xgboost","gbm")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(glmnet)
  library(Matrix)
  library(matrixStats)
  library(e1071)
  library(ranger)
  library(xgboost)
  library(gbm)
})

# ---------------------------
# Paths
# ---------------------------
DIR_MOD <- if (exists("DIR_MOD")) DIR_MOD else file.path("results","models")
dir.create(DIR_MOD, showWarnings=FALSE, recursive=TRUE)

f_metrics   <- file.path(DIR_MOD, "STRICT_LODO_multimodel_metrics.tsv")
f_pred_long <- file.path(DIR_MOD, "STRICT_LODO_multimodel_predictions_long.tsv")
f_imp_long  <- file.path(DIR_MOD, "STRICT_LODO_multimodel_importance_long.tsv")
f_rank      <- file.path(DIR_MOD, "STRICT_LODO_gene_ranking_aggregated.tsv")

# ---------------------------
# Helpers: IO + checks
# ---------------------------
must_exist <- function(f) if (!file.exists(f)) stop("Missing file: ", f)

read_limma_table <- function(ds) {
  f <- file.path(DIR_TAB, paste0(ds, "_limma_all.tsv"))
  must_exist(f)
  tt <- fread(f)
  
  gene_col  <- c("gene","Gene","SYMBOL","symbol")[c("gene","Gene","SYMBOL","symbol") %in% names(tt)][1]
  if (is.na(gene_col)) stop("No gene column in limma table: ", f)
  
  logfc_col <- c("logFC","log2FC","logFC_mrna")[c("logFC","log2FC","logFC_mrna") %in% names(tt)][1]
  if (is.na(logfc_col)) stop("No logFC column in limma table: ", f)
  
  fdr_candidates <- c("FDR","adj.P.Val","padj","adj_p","fdr","adjP")
  fdr_col <- fdr_candidates[fdr_candidates %in% names(tt)][1]
  if (is.na(fdr_col)) stop("No FDR/adj.P.Val column in limma table: ", f)
  
  out <- tt[, .(
    gene  = as.character(get(gene_col)),
    logFC = as.numeric(get(logfc_col)),
    fdr   = as.numeric(get(fdr_col))
  )]
  out <- out[is.finite(logFC) & is.finite(fdr)]
  out
}

load_expr_and_labels <- function(ds) {
  expr_path  <- file.path(DIR_PROC, paste0(ds, "_expr_gene.rds"))
  sheet_path <- file.path(DIR_META, paste0(ds, "_sample_sheet.csv"))
  must_exist(expr_path); must_exist(sheet_path)
  
  expr  <- readRDS(expr_path)   # genes x samples
  sheet <- fread(sheet_path)
  
  if (!("final_group" %in% names(sheet))) stop("final_group missing in ", sheet_path)
  if (!("sample_id" %in% names(sheet)))   stop("sample_id missing in ", sheet_path)
  
  sheet[final_group=="", final_group := NA_character_]
  sheet <- sheet[!is.na(final_group) & final_group %in% c("HGSC","Normal")]
  
  sids <- intersect(colnames(expr), sheet$sample_id)
  if (length(sids) < 6) stop("Too few matched samples in ", ds)
  
  sheet <- sheet[match(sids, sample_id)]
  expr  <- expr[, sids, drop=FALSE]
  
  y <- ifelse(sheet$final_group=="HGSC", 1L, 0L)
  list(expr=expr, y=y, sample_id=sids)
}

# ---------------------------
# TRAIN-only candidate pool (Top-N)
# strict: all train datasets must pass fdr_cut AND same direction.
# score: (-log10(mean_fdr)) * mean_abs_logFC
# ---------------------------
candidate_pool_train_only_topN <- function(train_ds, top_n=500, fdr_cut=0.05, min_abs_logFC=0) {
  tabs <- lapply(train_ds, function(ds) {
    tt <- read_limma_table(ds)
    setnames(tt, c("logFC","fdr"), c(paste0("logFC_", ds), paste0("fdr_", ds)))
    tt
  })
  
  common <- Reduce(intersect, lapply(tabs, `[[`, "gene"))
  if (length(common) < 50) return(character(0))
  
  tabs <- lapply(tabs, function(x) x[gene %in% common])
  m <- Reduce(function(a,b) merge(a,b, by="gene", all=FALSE), tabs)
  
  logfc_cols <- paste0("logFC_", train_ds)
  fdr_cols   <- paste0("fdr_",   train_ds)
  
  keep <- rep(TRUE, nrow(m))
  for (j in seq_along(train_ds)) {
    keep <- keep & (m[[fdr_cols[j]]] < fdr_cut) & (abs(m[[logfc_cols[j]]]) >= min_abs_logFC)
  }
  if (!any(keep)) return(character(0))
  mm <- m[keep]
  
  dirs <- sign(as.matrix(mm[, ..logfc_cols]))
  same_dir <- apply(dirs, 1, function(v) length(unique(v))==1)
  mm <- mm[same_dir]
  if (nrow(mm) == 0) return(character(0))
  
  mm[, mean_fdr := rowMeans(as.matrix(.SD)), .SDcols=fdr_cols]
  mm[, mean_abs_lfc := rowMeans(abs(as.matrix(.SD))), .SDcols=logfc_cols]
  mm[, score := (-log10(mean_fdr + 1e-12)) * mean_abs_lfc]
  
  setorder(mm, -score)
  head(mm$gene, min(top_n, nrow(mm)))
}

# ---------------------------
# Train-only scaling
# ---------------------------
scale_train_apply_test <- function(Xtr, Xte) {
  mu <- colMeans(Xtr)
  sd <- apply(Xtr, 2, sd)
  sd[!is.finite(sd) | sd==0] <- 1
  
  Xtr_s <- sweep(sweep(Xtr, 2, mu, "-"), 2, sd, "/")
  Xte_s <- sweep(sweep(Xte, 2, mu, "-"), 2, sd, "/")
  list(Xtr=Xtr_s, Xte=Xte_s)
}

# ---------------------------
# Metrics
# ---------------------------
auc_metrics <- function(y, prob, thr=0.5) {
  r <- pROC::roc(y, prob, quiet=TRUE, direction="<")
  auc <- as.numeric(pROC::auc(r))
  
  pred <- as.integer(prob >= thr)
  TP <- sum(y==1 & pred==1); TN <- sum(y==0 & pred==0)
  FP <- sum(y==0 & pred==1); FN <- sum(y==1 & pred==0)
  acc <- (TP+TN)/length(y)
  sens <- ifelse((TP+FN)==0, NA, TP/(TP+FN))
  spec <- ifelse((TN+FP)==0, NA, TN/(TN+FP))
  prec <- ifelse((TP+FP)==0, NA, TP/(TP+FP))
  f1 <- ifelse(is.na(prec) | is.na(sens) | (prec+sens)==0, NA, 2*prec*sens/(prec+sens))
  
  data.table(AUC=auc, TP=TP, TN=TN, FP=FP, FN=FN, Accuracy=acc, Sensitivity=sens, Specificity=spec, Precision=prec, F1=f1)
}

make_strat_foldid <- function(y, k) {
  y <- as.integer(y)
  idx1 <- which(y==1); idx0 <- which(y==0)
  if (length(idx1) < 2 || length(idx0) < 2) return(NULL)
  
  foldid <- integer(length(y))
  idx1 <- sample(idx1); idx0 <- sample(idx0)
  foldid[idx1] <- rep(1:k, length.out=length(idx1))
  foldid[idx0] <- rep(1:k, length.out=length(idx0))
  foldid
}

# ---------------------------
# Model trainers + importance
# IMPORTANT: All models will see SAFE gene names;
# we will map importance back to ORIGINAL gene symbols afterwards.
# ---------------------------
fit_elasticnet <- function(Xtr, ytr, alpha=0.5) {
  ytr <- as.integer(ytr)
  n0 <- sum(ytr==0); n1 <- sum(ytr==1)
  k <- min(5, floor(min(n0,n1)))
  if (!is.finite(k) || k < 2) k <- 2
  foldid <- make_strat_foldid(ytr, k)
  tm <- if (min(n0,n1) >= 10 && k >= 3) "auc" else "deviance"
  
  cv <- cv.glmnet(
    x=Xtr, y=ytr, family="binomial", alpha=alpha,
    nfolds=k, foldid=foldid, type.measure=tm
  )
  fit <- glmnet(Xtr, ytr, family="binomial", alpha=alpha, lambda=cv$lambda.1se)
  
  beta <- as.matrix(coef(fit))[-1,1]
  imp <- abs(beta)
  names(imp) <- colnames(Xtr)
  list(model=fit, importance=imp)
}
pred_elasticnet <- function(obj, Xte) as.numeric(predict(obj$model, Xte, type="response"))

fit_svm_linear <- function(Xtr, ytr) {
  yfac <- factor(ytr, levels=c(0,1))
  fit <- e1071::svm(
    x=Xtr, y=yfac,
    kernel="linear", probability=TRUE, scale=FALSE,
    type="C-classification"
  )
  w <- tryCatch({
    as.numeric(t(fit$coefs) %*% fit$SV)
  }, error=function(e) rep(0, ncol(Xtr)))
  names(w) <- colnames(Xtr)
  imp <- abs(w)
  list(model=fit, importance=imp)
}
pred_svm_linear <- function(obj, Xte) {
  pr <- predict(obj$model, Xte, probability=TRUE)
  probs <- attr(pr, "probabilities")
  if ("1" %in% colnames(probs)) return(as.numeric(probs[, "1"]))
  as.numeric(probs[, ncol(probs)])
}

fit_rf <- function(Xtr, ytr) {
  # Now safe colnames => formula interface OK
  df <- data.frame(y=factor(ytr, levels=c(0,1)), Xtr, check.names=FALSE)
  fit <- ranger(
    y ~ ., data=df,
    probability=TRUE,
    num.trees=1000,
    mtry=max(1, floor(sqrt(ncol(Xtr)))),
    importance="impurity",
    seed=1
  )
  imp <- fit$variable.importance
  imp <- imp[colnames(Xtr)]
  imp[is.na(imp)] <- 0
  list(model=fit, importance=as.numeric(imp))
}
pred_rf <- function(obj, Xte) {
  p <- predict(obj$model, data=data.frame(Xte, check.names=FALSE))$predictions
  as.numeric(p[, "1"])
}

fit_xgb <- function(Xtr, ytr) {
  ytr <- as.numeric(ytr)
  n0 <- sum(ytr==0); n1 <- sum(ytr==1)
  
  dtrain <- xgb.DMatrix(data=as.matrix(Xtr), label=ytr)
  
  # If data is scarce, AUC might not be computable in CV => logloss is safer
  use_auc <- (min(n0, n1) >= 10)
  eval_metric <- if (use_auc) "auc" else "logloss"
  
  nfold <- min(5, floor(min(n0, n1)))
  if (!is.finite(nfold) || nfold < 2) nfold <- 2
  
  params <- list(
    objective="binary:logistic",
    eval_metric=eval_metric,
    max_depth=3,
    eta=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    lambda=1
  )
  
  cv <- xgb.cv(
    params=params,
    data=dtrain,
    nrounds=300,
    nfold=nfold,
    stratified=TRUE,
    verbose=0,
    early_stopping_rounds=25
  )
  
  # best_n robust
  best_n <- cv$best_iteration
  if (length(best_n) == 0 || !is.finite(best_n)) {
    el <- cv$evaluation_log
    if (!is.null(el) && nrow(el) > 0) {
      if (eval_metric == "auc") {
        col <- grep("test.*auc.*mean", names(el), value=TRUE)
        if (length(col) >= 1) best_n <- which.max(el[[col[1]]])
      } else {
        col <- grep("test.*logloss.*mean", names(el), value=TRUE)
        if (length(col) >= 1) best_n <- which.min(el[[col[1]]])
      }
      if (!is.finite(best_n) || length(best_n)==0) best_n <- nrow(el)
    } else {
      best_n <- 200L
    }
  }
  
  best_n <- as.integer(best_n)
  if (!is.finite(best_n) || best_n < 10) best_n <- 200L
  
  fit <- xgb.train(params=params, data=dtrain, nrounds=best_n, verbose=0)
  
  imp_dt <- tryCatch(xgb.importance(model=fit), error=function(e) NULL)
  imp <- setNames(rep(0, ncol(Xtr)), colnames(Xtr))
  if (!is.null(imp_dt) && nrow(imp_dt) > 0) {
    imp[imp_dt$Feature] <- imp_dt$Gain
  }
  
  list(model=fit, best_n=best_n, importance=imp)
}

pred_xgb <- function(obj, Xte) {
  as.numeric(predict(obj$model, xgb.DMatrix(as.matrix(Xte))))
}
# ---------------------------
# GBM (FIX MISSING FUNCTION)
# ---------------------------
fit_gbm <- function(Xtr, ytr) {
  df <- data.frame(y=as.numeric(ytr), Xtr, check.names=FALSE)
  
  fit <- gbm::gbm(
    y ~ ., data=df,
    distribution="bernoulli",
    n.trees=200,
    interaction.depth=3,
    shrinkage=0.05,
    n.minobsinnode=10,
    verbose=FALSE
  )
  
  imp_df <- summary(fit, plotit=FALSE)
  imp <- setNames(rep(0, ncol(Xtr)), colnames(Xtr))
  
  if (!is.null(imp_df) && nrow(imp_df) > 0) {
    imp[imp_df$var] <- imp_df$rel.inf
  }
  
  list(model=fit, importance=imp)
}

pred_gbm <- function(obj, Xte) {
  as.numeric(
    predict(
      obj$model,
      newdata=data.frame(Xte, check.names=FALSE),
      n.trees=obj$model$n.trees,
      type="response"
    )
  )
}

# ---------------------------
# Datasets
# ---------------------------
datasets <- GSE_mrna
datasets <- datasets[
  file.exists(file.path(DIR_PROC, paste0(datasets, "_expr_gene.rds"))) &
    file.exists(file.path(DIR_META, paste0(datasets, "_sample_sheet.csv"))) &
    file.exists(file.path(DIR_TAB,  paste0(datasets, "_limma_all.tsv")))
]
if (length(datasets) < 3) stop("Need >=3 datasets with expr+sheet+limma_all.tsv")

DL <- lapply(datasets, load_expr_and_labels)
names(DL) <- datasets

# ---------------------------
# Run LODO
# ---------------------------
set.seed(1)

models <- c("ElasticNet","SVM","RandomForest","XGBoost","GBM")

all_metrics <- list()
all_preds   <- list()
all_imps    <- list()

TOP_N <- 500

for (test_ds in datasets) {
  
  train_ds <- setdiff(datasets, test_ds)
  
  cand <- candidate_pool_train_only_topN(train_ds, top_n=TOP_N, fdr_cut=0.05)
  if (length(cand) < 100) {
    message("[INFO] Few strict candidates in test=", test_ds, " (n=", length(cand), "). Trying relaxed fdr<0.10")
    cand2 <- candidate_pool_train_only_topN(train_ds, top_n=TOP_N, fdr_cut=0.10)
    if (length(cand2) > length(cand)) cand <- cand2
  }
  if (length(cand) < 50) stop("Candidate pool too small in test=", test_ds, " (n=", length(cand), ")")
  
  common_all <- Reduce(intersect, lapply(c(train_ds, test_ds), function(ds) rownames(DL[[ds]]$expr)))
  genes_use <- intersect(common_all, cand)
  if (length(genes_use) < 50) stop("genes_use too small after intersect in test=", test_ds, " (n=", length(genes_use), ")")
  
  # Build matrices: samples x genes
  Xtr0 <- do.call(rbind, lapply(train_ds, function(ds) t(DL[[ds]]$expr[genes_use, , drop=FALSE])))
  ytr  <- unlist(lapply(train_ds, function(ds) DL[[ds]]$y))
  Xte0 <- t(DL[[test_ds]]$expr[genes_use, , drop=FALSE])
  yte  <- DL[[test_ds]]$y
  
  # train-only scaling
  sc <- scale_train_apply_test(Xtr0, Xte0)
  Xtr <- sc$Xtr
  Xte <- sc$Xte
  
  # --------- CRITICAL FIX ----------
  # Make gene names SAFE for formula-based models (ranger/gbm),
  # but keep a map to ORIGINAL names for importance/ranking outputs.
  safe_names <- make.names(genes_use, unique=TRUE)
  safe_to_orig <- setNames(genes_use, safe_names)     # SAFE -> ORIGINAL
  orig_to_safe <- setNames(safe_names, genes_use)     # ORIGINAL -> SAFE
  
  colnames(Xtr) <- safe_names
  colnames(Xte) <- safe_names
  
  # ---- Fit models (TRAIN) and test on held-out dataset
  fit_list <- list()
  fit_list$ElasticNet    <- fit_elasticnet(Xtr, ytr, alpha=0.5)
  fit_list$SVM           <- fit_svm_linear(Xtr, ytr)
  fit_list$RandomForest  <- fit_rf(Xtr, ytr)
  fit_list$XGBoost       <- fit_xgb(Xtr, ytr)
  fit_list$GBM           <- fit_gbm(Xtr, ytr)
  
  pred_fun <- list(
    ElasticNet = pred_elasticnet,
    SVM        = pred_svm_linear,
    RandomForest = pred_rf,
    XGBoost      = pred_xgb,
    GBM          = pred_gbm
  )
  
  for (m in models) {
    prob <- pred_fun[[m]](fit_list[[m]], Xte)
    prob <- pmin(pmax(prob, 0), 1)
    
    met <- auc_metrics(yte, prob, thr=0.5)
    met[, `:=`(test_dataset=test_ds, train_datasets=paste(train_ds, collapse="|"), model=m, n_genes=length(genes_use), topN=TOP_N)]
    all_metrics[[length(all_metrics)+1]] <- met
    
    # predictions
    all_preds[[length(all_preds)+1]] <- data.table(
      test_dataset=test_ds,
      model=m,
      sample_id=DL[[test_ds]]$sample_id,
      y=yte,
      prob=prob
    )
    
    # importance (convert SAFE names -> ORIGINAL gene symbols)
    imp_safe <- fit_list[[m]]$importance
    if (is.null(names(imp_safe))) names(imp_safe) <- colnames(Xtr)
    imp_safe <- imp_safe[safe_names]
    imp_safe[is.na(imp_safe)] <- 0
    
    gene_orig <- safe_to_orig[names(imp_safe)]
    all_imps[[length(all_imps)+1]] <- data.table(
      test_dataset=test_ds,
      model=m,
      gene=as.character(gene_orig),
      importance=as.numeric(imp_safe)
    )
    
    message(sprintf("Fold test=%s | %-12s | AUC=%.3f | genes=%d", test_ds, m, met$AUC, length(genes_use)))
  }
}

metrics_dt <- rbindlist(all_metrics, fill=TRUE)
pred_dt    <- rbindlist(all_preds, fill=TRUE)
imp_dt     <- rbindlist(all_imps, fill=TRUE)

fwrite(metrics_dt, f_metrics, sep="\t")
fwrite(pred_dt,    f_pred_long, sep="\t")
fwrite(imp_dt,     f_imp_long, sep="\t")

message("\nSaved:")
message(" - ", f_metrics)
message(" - ", f_pred_long)
message(" - ", f_imp_long)

# ---------------------------
# Aggregate gene ranking across folds+models
# Reciprocal-rank score: mean(1/rank)
# ---------------------------
imp_dt[, rank := frank(-importance, ties.method="average"), by=.(test_dataset, model)]
imp_dt[, rr := 1 / rank]

rank_dt <- imp_dt[, .(
  mean_rr = mean(rr, na.rm=TRUE),
  median_rank = median(rank, na.rm=TRUE),
  mean_importance = mean(importance, na.rm=TRUE),
  n_entries = .N
), by=.(gene)]

setorder(rank_dt, -mean_rr, median_rank)

# OPTIONAL: annotate with global META if exists (annotation only; not used in training)
meta_f <- file.path(DIR_TAB, "META_consensus_STRICT.tsv")
if (!file.exists(meta_f)) meta_f <- file.path(DIR_TAB, "META_consensus_RELAXED.tsv")

if (file.exists(meta_f)) {
  meta <- fread(meta_f)
  gene_col <- c("gene","Gene","SYMBOL","symbol")[c("gene","Gene","SYMBOL","symbol") %in% names(meta)][1]
  logfc_col <- c("meta_logFC","meta_log2FC","meta_logfc","meta_LFC")[c("meta_logFC","meta_log2FC","meta_logfc","meta_LFC") %in% names(meta)][1]
  fdr_col <- c("meta_fdr","meta_FDR","FDR","padj","adj.P.Val","adj_p")[c("meta_fdr","meta_FDR","FDR","padj","adj.P.Val","adj_p") %in% names(meta)][1]
  
  if (!is.na(gene_col)) meta[, gene := as.character(get(gene_col))]
  if (!is.na(logfc_col)) meta[, meta_logFC := as.numeric(get(logfc_col))]
  if (!is.na(fdr_col)) meta[, meta_fdr := as.numeric(get(fdr_col))]
  
  meta2 <- unique(meta[, .(gene, meta_logFC, meta_fdr)])
  rank_dt <- merge(rank_dt, meta2, by="gene", all.x=TRUE)
  rank_dt[, direction := fifelse(is.finite(meta_logFC) & meta_logFC>0, "UP",
                                 fifelse(is.finite(meta_logFC) & meta_logFC<0, "DOWN", NA_character_))]
}

fwrite(rank_dt, f_rank, sep="\t")
message(" - ", f_rank)

message("\nTop 20 ranked genes:")
print(head(rank_dt, 20))

