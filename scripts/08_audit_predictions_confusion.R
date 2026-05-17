# scripts/08_audit_predictions_confusion.R
source("scripts/01_params.R")

pkgs <- c("data.table","pROC")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
suppressPackageStartupMessages({library(data.table); library(pROC)})

pred_file <- file.path(DIR_MOD, "LODO_nested_predictions.tsv")
stopifnot(file.exists(pred_file))
pred <- fread(pred_file)

# y باید 0/1 باشد
pred[, y := as.integer(y)]
pred[, pred_05 := ifelse(prob >= 0.5, 1L, 0L)]

# ---- sanity checks (label swap) ----
cat("\n[Check 1] Class counts per dataset:\n")
print(pred[, .N, by=.(test_dataset, y)])

cat("\n[Check 2] Mean prob by true class (should be higher for y=1):\n")
print(pred[, .(mean_prob_y1 = mean(prob[y==1], na.rm=TRUE),
               mean_prob_y0 = mean(prob[y==0], na.rm=TRUE)), by=test_dataset])

# ---- confusion matrix helper ----
conf_one <- function(dd) {
  TP <- sum(dd$y==1 & dd$pred_05==1)
  TN <- sum(dd$y==0 & dd$pred_05==0)
  FP <- sum(dd$y==0 & dd$pred_05==1)
  FN <- sum(dd$y==1 & dd$pred_05==0)
  acc <- (TP+TN)/max(1,(TP+TN+FP+FN))
  sens <- TP/max(1,(TP+FN))
  spec <- TN/max(1,(TN+FP))
  prec <- TP/max(1,(TP+FP))
  f1 <- if ((prec+sens)==0) NA else 2*prec*sens/(prec+sens)
  data.frame(TP=TP,TN=TN,FP=FP,FN=FN,Accuracy=acc,Sensitivity=sens,Specificity=spec,Precision=prec,F1=f1)
}

# ---- per dataset AUC + confusion ----
rows <- list()
for (ds in sort(unique(pred$test_dataset))) {
  dd <- pred[test_dataset==ds]
  auc <- as.numeric(pROC::auc(pROC::roc(dd$y, dd$prob, quiet=TRUE)))
  cf <- conf_one(dd)
  rows[[ds]] <- cbind(data.frame(test_dataset=ds, AUC=auc), cf)
}
out <- rbindlist(rows, fill=TRUE)
print(out)

dir.create(DIR_MOD, showWarnings = FALSE, recursive = TRUE)
f_out <- file.path(DIR_MOD, "Audit_confusion_by_dataset_threshold0.5.csv")
fwrite(out, f_out)
cat("\nSaved:", f_out, "\n")
