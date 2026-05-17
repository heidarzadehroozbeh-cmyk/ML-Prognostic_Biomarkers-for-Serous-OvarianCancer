# scripts/10_benchmark_models_STRICT_LODO.R
# Benchmark models on STRICT LODO panels (train-only scaling; no leakage)
# Models: ElasticNet (glmnet), RandomForest (ranger), XGBoost (xgboost)

source("scripts/01_params.R")

pkgs <- c("data.table","pROC","glmnet","Matrix","ranger","xgboost")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(glmnet)
  library(Matrix)
  library(ranger)
  library(xgboost)
})

# ---------- paths ----------
DIR_FIG <- if (exists("DIR_FIG")) DIR_FIG else "figures"
DIR_MOD <- if (exists("DIR_MOD")) DIR_MOD else file.path("results","models")
dir.create(DIR_MOD, showWarnings=FALSE, recursive=TRUE)

pred_f  <- file.path(DIR_MOD, "LODO_STRICT_predictions.tsv")
panel_f <- file.path(DIR_MOD, "LODO_STRICT_panels.tsv")

if (!file.exists(pred_f)) stop("Missing: ", pred_f)
if (!file.exists(panel_f)) stop("Missing: ", panel_f)

pred0  <- fread(pred_f)
panels <- fread(panel_f)

datasets <- sort(unique(pred0$test_dataset))
if (length(datasets) < 2) stop("Need >=2 datasets in predictions file.")

# ---------- helpers ----------
load_expr_and_sheet <- function(ds) {
  expr_path  <- file.path(DIR_PROC, paste0(ds, "_expr_gene.rds"))
  sheet_path <- file.path(DIR_META, paste0(ds, "_sample_sheet.csv"))
  if (!file.exists(expr_path)) stop("Missing expr: ", expr_path)
  if (!file.exists(sheet_path)) stop("Missing sheet: ", sheet_path)
  
  expr  <- readRDS(expr_path)   # genes x samples
  sheet <- fread(sheet_path)
  
  # expect columns: sample_id + final_group
  if (!("sample_id" %in% names(sheet))) stop("sample_id missing in ", sheet_path)
  if (!("final_group" %in% names(sheet))) stop("final_group missing in ", sheet_path)
  
  sheet$final_group[sheet$final_group==""] <- NA
  sheet <- sheet[!is.na(final_group) & final_group %in% c("HGSC","Normal")]
  
  sids <- intersect(colnames(expr), sheet$sample_id)
  if (length(sids) < 4) stop("Too few matched samples in ", ds, " (", length(sids), ")")
  
  sheet <- sheet[match(sids, sheet$sample_id)]
  expr  <- expr[, sids, drop=FALSE]
  
  list(expr=expr, sheet=sheet)
}

train_test_mats <- function(train_ds, test_ds, genes) {
  # combine TRAIN datasets
  Xtr_list <- list()
  ytr_list <- list()
  
  for (ds in train_ds) {
    o <- load_expr_and_sheet(ds)
    expr <- o$expr
    ss   <- o$sheet
    g2 <- intersect(genes, rownames(expr))
    if (length(g2) < 5) next
    X <- t(expr[g2, , drop=FALSE])  # samples x genes
    y <- ifelse(ss$final_group=="HGSC", 1L, 0L)
    Xtr_list[[ds]] <- X
    ytr_list[[ds]] <- y
  }
  
  if (length(Xtr_list) == 0) stop("No train matrices built. Check genes overlap.")
  Xtr <- do.call(rbind, Xtr_list)
  ytr <- as.integer(unlist(ytr_list))
  
  # TEST dataset
  ot <- load_expr_and_sheet(test_ds)
  exprt <- ot$expr
  sst   <- ot$sheet
  g2t <- intersect(genes, rownames(exprt))
  if (length(g2t) < 5) stop("Too few genes in test overlap: ", test_ds)
  
  Xte <- t(exprt[g2t, , drop=FALSE])
  yte <- as.integer(ifelse(sst$final_group=="HGSC", 1L, 0L))
  
  # align columns
  common <- intersect(colnames(Xtr), colnames(Xte))
  if (length(common) < 5) stop("Too few common genes after align: ", length(common))
  Xtr <- Xtr[, common, drop=FALSE]
  Xte <- Xte[, common, drop=FALSE]
  
  # ---- train-only scaling ----
  mu  <- colMeans(Xtr, na.rm=TRUE)
  sdv <- apply(Xtr, 2, sd, na.rm=TRUE)
  sdv[!is.finite(sdv) | sdv==0] <- 1
  
  Xtrz <- sweep(sweep(Xtr, 2, mu, "-"), 2, sdv, "/")
  Xtez <- sweep(sweep(Xte, 2, mu, "-"), 2, sdv, "/")
  
  list(Xtr=Xtrz, ytr=ytr, Xte=Xtez, yte=yte, genes=common)
}

select_panel_for_fold <- function(test_ds, freq_cut=0.60, min_genes=12, max_genes=80) {
  pp <- panels[test_dataset==test_ds]
  if (nrow(pp) == 0) stop("No rows for test_dataset in panels: ", test_ds)
  if (!("gene" %in% names(pp))) stop("panels file missing 'gene'")
  if (!("freq" %in% names(pp))) stop("panels file missing 'freq'")
  
  setorder(pp, -freq)
  
  g <- pp[freq >= freq_cut, gene]
  if (length(g) < min_genes) g <- pp[1:min(min_genes, .N), gene]
  if (length(g) > max_genes) g <- g[1:max_genes]
  unique(as.character(g))
}

auc_metrics <- function(y, prob, thr=0.5) {
  auc <- as.numeric(pROC::auc(pROC::roc(y, prob, quiet=TRUE, direction="<")))
  pred <- as.integer(prob >= thr)
  TP <- sum(y==1 & pred==1); TN <- sum(y==0 & pred==0)
  FP <- sum(y==0 & pred==1); FN <- sum(y==1 & pred==0)
  acc <- (TP+TN)/(TP+TN+FP+FN)
  sens <- ifelse((TP+FN)==0, NA, TP/(TP+FN))
  spec <- ifelse((TN+FP)==0, NA, TN/(TN+FP))
  brier <- mean((prob - y)^2)
  data.frame(AUC=auc, TP=TP,TN=TN,FP=FP,FN=FN, Accuracy=acc, Sensitivity=sens, Specificity=spec, Brier=brier)
}

# ---------- model trainers ----------
# UPDATED: dynamic nfolds + stratified foldid to keep type.measure="auc" when possible
fit_glmnet <- function(X, y) {
  x <- as.matrix(X)
  y <- as.integer(y)
  n <- nrow(x)
  
  # want >=10 obs per fold for AUC
  nfolds_auc <- floor(n / 10)
  
  if (nfolds_auc >= 3) {
    nfolds <- min(5, nfolds_auc)  # 3..5
    measure <- "auc"
  } else {
    nfolds <- 3
    measure <- "deviance"
  }
  
  set.seed(2025)
  foldid <- integer(n)
  for (cls in c(0L, 1L)) {
    idx <- which(y == cls)
    foldid[idx] <- sample(rep(seq_len(nfolds), length.out = length(idx)))
  }
  
  cv.glmnet(
    x, y,
    family="binomial",
    alpha=0.5,
    type.measure=measure,
    standardize=FALSE,
    nfolds=nfolds,
    foldid=foldid
  )
}

pred_glmnet <- function(fit, X) as.numeric(predict(fit, as.matrix(X), s="lambda.1se", type="response"))

fit_rf <- function(X, y) {
  df <- data.frame(y=factor(y, levels=c(0,1)), X, check.names=FALSE)
  ranger(y ~ ., data=df, probability=TRUE, num.trees=800,
         mtry=max(1, floor(sqrt(ncol(X)))))
}
pred_rf <- function(fit, X) {
  p <- predict(fit, data=data.frame(X, check.names=FALSE))$predictions
  as.numeric(p[, "1"])
}

fit_xgb <- function(X, y) {
  dtrain <- xgb.DMatrix(data=as.matrix(X), label=as.numeric(y))
  params <- list(
    objective="binary:logistic",
    eval_metric="auc",
    max_depth=3,
    eta=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    min_child_weight=1,
    lambda=1
  )
  xgb.train(params=params, data=dtrain, nrounds=250, verbose=0)
}
pred_xgb <- function(fit, X) as.numeric(predict(fit, xgb.DMatrix(as.matrix(X))))

# ---------- run LODO benchmark ----------
all_pred <- list()
all_sum  <- list()

models <- c("ElasticNet","RandomForest","XGBoost")

for (test_ds in datasets) {
  train_ds <- setdiff(datasets, test_ds)
  
  panel_genes <- select_panel_for_fold(test_ds, freq_cut=0.60, min_genes=12, max_genes=80)
  dat <- train_test_mats(train_ds, test_ds, panel_genes)
  
  # ElasticNet
  fit1 <- fit_glmnet(dat$Xtr, dat$ytr)
  pr1  <- pred_glmnet(fit1, dat$Xte)
  
  # RF
  fit2 <- fit_rf(dat$Xtr, dat$ytr)
  pr2  <- pred_rf(fit2, dat$Xte)
  
  # XGB
  fit3 <- fit_xgb(dat$Xtr, dat$ytr)
  pr3  <- pred_xgb(fit3, dat$Xte)
  
  dd <- data.frame(
    test_dataset=test_ds,
    y=dat$yte,
    ElasticNet=pr1,
    RandomForest=pr2,
    XGBoost=pr3,
    stringsAsFactors=FALSE
  )
  all_pred[[test_ds]] <- dd
  
  for (m in models) {
    met <- auc_metrics(dd$y, dd[[m]], thr=0.5)
    met$model <- m
    met$test_dataset <- test_ds
    met$panel_size <- ncol(dat$Xte)
    all_sum[[paste(test_ds,m,sep="__")]] <- met
  }
  
  message("Done fold: test=", test_ds, " | panel=", ncol(dat$Xte))
}

pred_out <- rbindlist(all_pred, fill=TRUE)
sum_out  <- rbindlist(all_sum, fill=TRUE)

f_pred <- file.path(DIR_MOD, "LODO_STRICT_benchmark_predictions_wide.tsv")
f_sum  <- file.path(DIR_MOD, "LODO_STRICT_benchmark_metrics.tsv")
fwrite(pred_out, f_pred, sep="\t")
fwrite(sum_out,  f_sum,  sep="\t")

message("Saved: ", f_pred)
message("Saved: ", f_sum)
