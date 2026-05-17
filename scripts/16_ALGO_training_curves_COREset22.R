# scripts/16_ALGO_training_curves_COREset22.R
# Training/Tuning curves (NOT ROC) for: GBM, SVM, ElasticNet, XGBoost, RandomForest
# Feature space: your 22 meta-analysis genes (from figures/FINAL22_gene_list_*.tsv)

source("scripts/01_params.R")

# ---------- packages ----------
pkgs <- c("data.table","ggplot2","patchwork","Matrix","pROC","glmnet")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(Matrix); library(pROC); library(glmnet)
})

safe_has <- function(pkg) requireNamespace(pkg, quietly=TRUE)

DIR_FIG <- if (exists("DIR_FIG")) DIR_FIG else "figures"
dir.create(DIR_FIG, showWarnings=FALSE, recursive=TRUE)

must_exist <- function(f) if (!file.exists(f)) stop("Missing file: ", f)

# ---------- inputs ----------
# 22 meta genes list (your file already exists according to your console)
final22_files <- list.files("figures", pattern="FINAL22_gene_list_.*\\.tsv$", full.names=TRUE)
if (length(final22_files) == 0) stop("Cannot find FINAL22_gene_list_*.tsv in ./figures/")
FINAL22_F <- final22_files[1]
message("[INFO] Using gene list: ", FINAL22_F)

gdt <- fread(FINAL22_F)
if (!("gene" %in% names(gdt))) stop("FINAL22 file must contain column: gene")
GENES <- unique(gdt$gene)
if (length(GENES) < 5) stop("Too few genes in FINAL22 list.")

# datasets detected from DEG tables
deg_files <- list.files("results/tables", pattern="_limma_DEG\\.tsv$", full.names=TRUE)
if (length(deg_files) == 0) stop("No DEG files found: results/tables/*_limma_DEG.tsv")
DSETS0 <- sort(unique(gsub("_limma_DEG\\.tsv$", "", basename(deg_files))))

has_inputs <- function(ds) {
  file.exists(file.path(DIR_PROC, paste0(ds, "_expr_gene.rds"))) &&
    file.exists(file.path(DIR_META, paste0(ds, "_sample_sheet.csv")))
}
DSETS <- DSETS0[vapply(DSETS0, has_inputs, logical(1))]
if (length(DSETS) < 2) stop("Need >=2 datasets with expr+sheet available. Found: ", paste(DSETS, collapse=", "))

message("[INFO] Datasets used: ", paste(DSETS, collapse=", "))

# ---------- load pooled X,y using within-dataset z-scoring ----------
load_xy <- function(ds, genes) {
  expr_path  <- file.path(DIR_PROC, paste0(ds, "_expr_gene.rds"))
  sheet_path <- file.path(DIR_META, paste0(ds, "_sample_sheet.csv"))
  must_exist(expr_path); must_exist(sheet_path)
  
  expr  <- readRDS(expr_path)   # genes x samples
  sheet <- fread(sheet_path)
  
  if (!all(c("sample_id","final_group") %in% names(sheet))) stop("sample_sheet missing sample_id/final_group: ", sheet_path)
  sheet[final_group=="", final_group := NA_character_]
  sheet <- sheet[final_group %in% c("HGSC","Normal")]
  sids <- intersect(colnames(expr), sheet$sample_id)
  sheet <- sheet[match(sids, sample_id)]
  expr  <- expr[, sids, drop=FALSE]
  
  g <- intersect(genes, rownames(expr))
  if (length(g) < 5) stop("Too few FINAL22 genes present in: ", ds, " (found ", length(g), ")")
  
  X <- t(expr[g, , drop=FALSE])  # samples x genes
  
  # within-dataset z-score per gene
  Xz <- apply(X, 2, function(v){
    m <- mean(v, na.rm=TRUE); s <- sd(v, na.rm=TRUE)
    if (!is.finite(s) || s==0) return(v*0)
    (v-m)/s
  })
  
  y <- ifelse(sheet$final_group=="HGSC", 1L, 0L)
  list(X=Xz, y=y, ds=ds, genes=colnames(Xz))
}

DL <- lapply(DSETS, load_xy, genes=GENES)
names(DL) <- DSETS

# use common genes across all datasets (safe)
common_genes <- Reduce(intersect, lapply(DL, function(z) colnames(z$X)))
if (length(common_genes) < 5) stop("Common genes across datasets < 5. Common genes: ", paste(common_genes, collapse=", "))

X_list <- lapply(DL, function(z) z$X[, common_genes, drop=FALSE])
y_list <- lapply(DL, function(z) z$y)
X <- do.call(rbind, X_list)
y <- unlist(y_list)

# basic sanity
message("[INFO] Pooled samples: ", nrow(X), " | Features: ", ncol(X))

set.seed(1)

# K-fold CV splits for methods where we implement CV manually
K <- 5
fold_id <- sample(rep(1:K, length.out=length(y)))

auc01 <- function(y_true, prob) {
  r <- pROC::roc(y_true, prob, quiet=TRUE, direction="<")
  as.numeric(pROC::auc(r))
}

# ---------- Panel A: ElasticNet CV curve (AUC vs log(lambda)) ----------
p_elnet <- {
  cvfit <- cv.glmnet(as.matrix(X), y, family="binomial", alpha=0.5,
                     nfolds=K, foldid=fold_id, type.measure="auc")
  cvdf <- data.table(lambda=cvfit$lambda, cvm=cvfit$cvm, cvs=cvfit$cvsd)
  ggplot(cvdf, aes(x=log(lambda), y=cvm)) +
    geom_line(linewidth=0.7) +
    geom_vline(xintercept=log(cvfit$lambda.1se), linetype="dashed", linewidth=0.5) +
    theme_bw(base_size=11) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(size=12, face="bold")) +
    labs(title="ElasticNet: CV curve", x="log(lambda)", y="CV AUC")
}

# ---------- Panel B: SVM (RBF) CV AUC vs log2(C) (gamma fixed = 1/p) ----------
p_svm <- {
  if (!safe_has("e1071")) {
    ggplot() + theme_void() + ggtitle("SVM: package e1071 not installed")
  } else {
    library(e1071)
    gamma0 <- 1 / ncol(X)
    C_grid <- 2^seq(-4, 8, by=1)
    
    res <- rbindlist(lapply(C_grid, function(Cv){
      aucs <- numeric(K)
      for (k in 1:K) {
        tr <- which(fold_id != k); te <- which(fold_id == k)
        fit <- e1071::svm(x=X[tr, , drop=FALSE], y=as.factor(y[tr]),
                          type="C-classification", kernel="radial",
                          cost=Cv, gamma=gamma0, probability=TRUE, scale=FALSE)
        pr <- attr(predict(fit, X[te, , drop=FALSE], probability=TRUE), "probabilities")
        # probability column name depends on factor levels; take "1" if exists else last column
        p1 <- if ("1" %in% colnames(pr)) pr[, "1"] else pr[, ncol(pr)]
        aucs[k] <- auc01(y[te], p1)
      }
      data.table(log2C=log2(Cv), AUC=mean(aucs, na.rm=TRUE))
    }))
    
    ggplot(res, aes(x=log2C, y=AUC)) +
      geom_line(linewidth=0.7) +
      geom_point(size=1.6) +
      theme_bw(base_size=11) +
      theme(panel.grid.minor=element_blank(),
            plot.title=element_text(size=12, face="bold")) +
      labs(title="SVM (RBF): CV curve", x="log2(C)", y="CV AUC")
  }
}

# ---------- Panel C: RandomForest OOB error vs #trees ----------
p_rf <- {
  if (!safe_has("randomForest")) {
    ggplot() + theme_void() + ggtitle("RandomForest: package randomForest not installed")
  } else {
    library(randomForest)
    # mtry default: sqrt(p)
    mtry0 <- max(1, floor(sqrt(ncol(X))))
    rf <- randomForest(x=X, y=as.factor(y), ntree=800, mtry=mtry0, importance=FALSE)
    err <- rf$err.rate[, "OOB"]
    dt <- data.table(trees=seq_along(err), OOB_error=as.numeric(err))
    ggplot(dt, aes(x=trees, y=OOB_error)) +
      geom_line(linewidth=0.7) +
      theme_bw(base_size=11) +
      theme(panel.grid.minor=element_blank(),
            plot.title=element_text(size=12, face="bold")) +
      labs(title=sprintf("RandomForest: OOB curve (mtry=%d)", mtry0),
           x="#Trees", y="OOB error")
  }
}

# ---------- Panel D: GBM CV deviance vs #trees ----------
p_gbm <- {
  if (!safe_has("gbm")) {
    ggplot() + theme_void() + ggtitle("GBM: package gbm not installed")
  } else {
    library(gbm)
    dt <- data.table(y=y, X)
    # gbm expects y numeric 0/1 for bernoulli
    fit <- gbm::gbm(
      y ~ ., data=dt,
      distribution="bernoulli",
      n.trees=2000,
      interaction.depth=3,
      shrinkage=0.01,
      n.minobsinnode=10,
      bag.fraction=0.7,
      cv.folds=K,
      verbose=FALSE
    )
    # cv.error length = n.trees
    cverr <- fit$cv.error
    d <- data.table(trees=seq_along(cverr), cv_deviance=as.numeric(cverr))
    ggplot(d, aes(x=trees, y=cv_deviance)) +
      geom_line(linewidth=0.7) +
      theme_bw(base_size=11) +
      theme(panel.grid.minor=element_blank(),
            plot.title=element_text(size=12, face="bold")) +
      labs(title="GBM: CV curve", x="#Trees", y="CV deviance (lower is better)")
  }
}

# ---------- Panel E: XGBoost CV AUC vs rounds ----------
p_xgb <- {
  if (!safe_has("xgboost")) {
    ggplot() + theme_void() + ggtitle("XGBoost: package xgboost not installed")
  } else {
    library(xgboost)
    dtrain <- xgb.DMatrix(data=as.matrix(X), label=y)
    params <- list(
      objective="binary:logistic",
      eval_metric="auc",
      eta=0.05,
      max_depth=4,
      subsample=0.8,
      colsample_bytree=0.8
    )
    cv <- xgb.cv(
      params=params,
      data=dtrain,
      nrounds=2000,
      nfold=K,
      stratified=TRUE,
      verbose=0,
      early_stopping_rounds=40
    )
    elog <- as.data.table(cv$evaluation_log)
    # columns typically: iter, train_auc_mean, test_auc_mean
    col_test <- grep("^test_auc_mean$", names(elog), value=TRUE)
    if (length(col_test)==0) col_test <- grep("test.*auc.*mean", names(elog), value=TRUE)[1]
    d <- data.table(round=elog$iter, test_auc_mean=elog[[col_test]])
    ggplot(d, aes(x=round, y=test_auc_mean)) +
      geom_line(linewidth=0.7) +
      theme_bw(base_size=11) +
      theme(panel.grid.minor=element_blank(),
            plot.title=element_text(size=12, face="bold")) +
      labs(title="XGBoost: CV curve", x="Boosting rounds", y="CV AUC")
  }
}

# ---------- assemble + save ----------
panel <- (p_gbm | p_svm) / (p_elnet | p_xgb) / p_rf +
  plot_annotation(
    title="Training / tuning curves for classical ML algorithms (22-gene meta feature space)",
    subtitle="Curves summarize internal CV/OOB behavior (diagnostics), not external test ROC",
    theme=theme(plot.title=element_text(size=14, face="bold"),
                plot.subtitle=element_text(size=11))
  )

OUT_PDF <- file.path(DIR_FIG, "ALGO_training_curves_panel.pdf")
OUT_PNG <- file.path(DIR_FIG, "ALGO_training_curves_panel.png")

ggsave(OUT_PDF, panel, width=12.0, height=14.0, device=cairo_pdf, limitsize=FALSE)
ggsave(OUT_PNG, panel, width=12.0, height=14.0, dpi=450, bg="white", limitsize=FALSE)

message("[OK] Saved: ", OUT_PDF)
message("[OK] Saved: ", OUT_PNG)

# Also save each panel separately (useful for PPT)
ggsave(file.path(DIR_FIG, "curve_ElasticNet_CV.pdf"), p_elnet, width=6.0, height=4.6, device=cairo_pdf, limitsize=FALSE)
ggsave(file.path(DIR_FIG, "curve_ElasticNet_CV.png"), p_elnet, width=6.0, height=4.6, dpi=600, bg="white", limitsize=FALSE)

ggsave(file.path(DIR_FIG, "curve_SVM_CV.pdf"), p_svm, width=6.0, height=4.6, device=cairo_pdf, limitsize=FALSE)
ggsave(file.path(DIR_FIG, "curve_SVM_CV.png"), p_svm, width=6.0, height=4.6, dpi=600, bg="white", limitsize=FALSE)

ggsave(file.path(DIR_FIG, "curve_RandomForest_OOB.pdf"), p_rf, width=6.0, height=4.6, device=cairo_pdf, limitsize=FALSE)
ggsave(file.path(DIR_FIG, "curve_RandomForest_OOB.png"), p_rf, width=6.0, height=4.6, dpi=600, bg="white", limitsize=FALSE)

ggsave(file.path(DIR_FIG, "curve_GBM_CV.pdf"), p_gbm, width=6.0, height=4.6, device=cairo_pdf, limitsize=FALSE)
ggsave(file.path(DIR_FIG, "curve_GBM_CV.png"), p_gbm, width=6.0, height=4.6, dpi=600, bg="white", limitsize=FALSE)

ggsave(file.path(DIR_FIG, "curve_XGBoost_CV.pdf"), p_xgb, width=6.0, height=4.6, device=cairo_pdf, limitsize=FALSE)
ggsave(file.path(DIR_FIG, "curve_XGBoost_CV.png"), p_xgb, width=6.0, height=4.6, dpi=600, bg="white", limitsize=FALSE)

message("DONE.")
