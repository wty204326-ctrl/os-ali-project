#!/usr/bin/env Rscript
# 05_validation.R — GSE111828 外部数据集验证诊断模型
# 数据: RNA-seq 逐样本 counts(xlsx combined 表)
# 跨平台:验证集基因独立 z-score(不套用芯片均值/标准差),再输入已在训练集 z-score 空间拟合的模型
# 用法: Rscript scripts/05_validation.R

ca <- commandArgs(trailingOnly = FALSE)
f  <- sub("^--file=", "", ca[grepl("^--file=", ca)])
script_dir <- if (length(f) == 1) dirname(normalizePath(f)) else "scripts"
source(file.path(script_dir, "_common.R"))

base_dir <- os_ali_base()
out_dir  <- file.path(base_dir, "results")
set.seed(20260909)
suppressMessages({
  library(readxl); library(glmnet); library(randomForest)
  library(xgboost); library(pROC); library(ggplot2)
})

## ---------- 1. 读取验证数据 ----------
f <- file.path(base_dir, "data/raw/GSE111828_results.xlsx")
d <- suppressMessages(read_excel(f, sheet = "combined"))
sample_cols <- grep("^MMR", names(d), value = TRUE)
stopifnot(length(sample_cols) == 24)
sym <- d$symbol
cnt <- as.matrix(d[, sample_cols])
rownames(cnt) <- ifelse(is.na(sym), paste0("NA_", seq_along(sym)), sym)
cnt <- cnt[!duplicated(rownames(cnt)), ]

## GEO: S01–S04 = 0h untreated control;其余为 12–72h APAP。按列名中的 S 编号分组,不依赖列顺序。
snum <- as.integer(sub("^.*_S([0-9]+)_.*$", "\\1", sample_cols))
if (any(is.na(snum))) {
  stop("无法从列名解析样本编号: ", paste(sample_cols[is.na(snum)], collapse = ", "))
}
group <- ifelse(snum <= 4, "Control", "APAP")
yv <- factor(group, levels = c("Control", "APAP"))
cat(sprintf("GSE111828: Control=%d  APAP=%d (按 S 编号,S01-S04=对照)\n",
            sum(group == "Control"), sum(group == "APAP")))

## ---------- 2. OS-DEGs 特征:log2 计数后按验证集独立 z-score ----------
mods <- readRDS(file.path(out_dir, "models_GSE17649.rds"))
feats <- as.character(mods$features)
missing <- setdiff(feats, rownames(cnt))
cat("特征基因在验证集中缺失:", ifelse(length(missing) == 0, "无", paste(missing, collapse = ",")), "\n")
present <- intersect(feats, rownames(cnt))
if (!length(present)) stop("验证集不含任何训练特征基因")
Xv <- t(log2(cnt[present, , drop = FALSE] + 1))
if (isTRUE(mods$scaled)) {
  Xv <- zscore_cols(Xv)
} else {
  cat("警告: 当前模型未在 z-score 空间训练。请先重新运行 scripts/03_ml_modeling.R。\n")
}
Xv <- align_features(Xv, feats)

if ("Hmox1" %in% colnames(Xv)) {
  hm <- Xv[, "Hmox1"]
  cat(sprintf("Hmox1 z-score 均值: Control=%.2f  APAP=%.2f (方向%s)\n",
              mean(hm[yv == "Control"]), mean(hm[yv == "APAP"]),
              ifelse(mean(hm[yv == "APAP"]) > mean(hm[yv == "Control"]), "正确", "异常")))
}

## ---------- 3. 三模型预测(列必须与训练特征一致) ----------
pv <- data.frame(sample = sample_cols, group = group, stringsAsFactors = FALSE)
pv$LASSO        <- as.numeric(predict(mods$lasso, newx = Xv,
                                      s = "lambda.min", type = "response"))
pv$randomForest <- predict(mods$rf, Xv, type = "prob")[, "APAP"]
pv$XGBoost      <- predict(mods$xgb, xgb.DMatrix(Xv))
write.csv(pv, file.path(out_dir, "validation_predictions.csv"), row.names = FALSE)

bin_metrics <- function(score, y, thr = 0.5) {
  pred <- factor(ifelse(score >= thr, "APAP", "Control"), levels = levels(y))
  tp <- sum(pred == "APAP" & y == "APAP")
  tn <- sum(pred == "Control" & y == "Control")
  fp <- sum(pred == "APAP" & y == "Control")
  fn <- sum(pred == "Control" & y == "APAP")
  list(sens = if ((tp + fn) == 0) NA_real_ else tp / (tp + fn),
       spec = if ((tn + fp) == 0) NA_real_ else tn / (tn + fp))
}

res <- data.frame()
roc_dat <- list()
for (m in c("LASSO", "randomForest", "XGBoost")) {
  rr <- roc(yv, pv[[m]], levels = c("Control", "APAP"), direction = "<", quiet = TRUE)
  bm <- bin_metrics(pv[[m]], yv, 0.5)
  res <- rbind(res, data.frame(model = m,
                               auc = round(as.numeric(auc(rr)), 4),
                               sens_0.5 = round(bm$sens, 4),
                               spec_0.5 = round(bm$spec, 4)))
  co <- coords(rr, "all", ret = c("specificity", "sensitivity"), as.list = FALSE)
  roc_dat[[m]] <- data.frame(model = m, fpr = 1 - co$specificity, tpr = co$sensitivity,
                             auc = round(as.numeric(auc(rr)), 3))
}
roc_all <- do.call(rbind, roc_dat)
write.csv(res, file.path(out_dir, "validation_auc.csv"), row.names = FALSE)
print(res, row.names = FALSE)

lab <- roc_all[, c("model", "auc")] |> unique()
p <- ggplot(roc_all, aes(fpr, tpr, color = model)) +
  geom_path(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_brewer(palette = "Set1",
                     labels = setNames(sprintf("%s (AUC=%.3f)", lab$model, lab$auc), lab$model)) +
  coord_equal() + theme_bw(base_size = 12) +
  labs(x = "1 − Specificity", y = "Sensitivity",
       title = "External validation ROC (GSE111828)",
       subtitle = "4 Control (0hr) vs 20 APAP (12-72hr); within-dataset z-score", color = NULL)
ggsave(file.path(out_dir, "roc_validation.png"), p, width = 6, height = 5.5, dpi = 300)
cat("[05] 外部验证完成\n")
