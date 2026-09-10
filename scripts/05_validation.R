#!/usr/bin/env Rscript
# 05_validation.R — GSE111828 外部数据集验证(区分度 + 校准检查)
#
# 方法学说明(修正 v1.0):
#   - 芯片 log2 强度与 RNA-seq log2 计数不在同一数值尺度,直接预测会导致
#     概率整体偏移(阈值不可迁移)。本版用训练集特征均值/标准差对验证集
#     做 z-score 标准化后再预测。
#   - AUC 衡量的是"区分度"(排序能力),不是校准的诊断概率;报告同时给出
#     两组预测概率的中位数与范围,阈值诊断需谨慎。
#   - 重复注释的基因取行均值(与训练集处理一致)。
#
# 用法: Rscript scripts/05_validation.R
# 输出: results/validation_predictions.csv
#       results/validation_auc.csv
#       results/validation_calibration.csv
#       results/roc_validation.png

sa <- commandArgs(FALSE)
fs <- sub("^--file=", "", sa[grepl("^--file=", sa)])
base_dir <- if (length(fs)) normalizePath(file.path(dirname(fs), "..")) else normalizePath(getwd())
out_dir <- file.path(base_dir, "results")
set.seed(20260909)
suppressMessages({
  library(readxl); library(glmnet); library(randomForest)
  library(xgboost); library(pROC); library(ggplot2)
})

## ---------- 1. 读取验证数据(24 样本逐基因 counts) ----------
f <- file.path(base_dir, "data/raw/GSE111828_results.xlsx")
d <- suppressMessages(read_excel(f, sheet = "combined"))
sample_cols <- grep("^MMR", names(d), value = TRUE)
stopifnot(length(sample_cols) == 24)
sym <- d$symbol
cnt <- as.matrix(d[, sample_cols])
rownames(cnt) <- ifelse(is.na(sym) | sym == "", paste0("NA_", seq_along(sym)), sym)

## 与训练集一致:重复基因取行均值
u <- unique(rownames(cnt))
cnt <- t(sapply(u, function(g) colMeans(cnt[rownames(cnt) == g, , drop = FALSE])))

group <- c(rep("Control", 4), rep("APAP", 20))
yv <- factor(group, levels = c("Control", "APAP"))

## ---------- 2. 特征与标准化(用训练集统计量;符号大小写不敏感) ----------
m <- readRDS(file.path(out_dir, "models_GSE17649.rds"))
feats <- toupper(m$features)
cnt_up <- cnt; rownames(cnt_up) <- toupper(rownames(cnt_up))
u <- unique(rownames(cnt_up))
cnt_up <- t(sapply(u, function(g) colMeans(cnt_up[rownames(cnt_up) == g, , drop = FALSE])))
missing <- setdiff(feats, rownames(cnt_up))
stopifnot(length(missing) == 0)   # 训练特征必须全部存在于验证数据
Xraw <- t(log2(cnt_up[feats, , drop = FALSE] + 1))
Zv <- sweep(sweep(Xraw, 2, m$train_mean), 2, m$train_sd, "/")

## Hmox1 方向合理性检查
if ("Hmox1" %in% feats) {
  hm <- Xraw[, "Hmox1"]
  cat(sprintf("Hmox1 log2cpm均值: Control=%.2f APAP=%.2f (方向%s)\n",
              mean(hm[yv == "Control"]), mean(hm[yv == "APAP"]),
              ifelse(mean(hm[yv == "APAP"]) > mean(hm[yv == "Control"]), "正确", "异常")))
}

## ---------- 3. 预测(标准化后) ----------
pv <- data.frame(sample = sample_cols, group = group)
pv$LASSO        <- as.numeric(predict(m$lasso, newx = Zv[, feats, drop = FALSE],
                                      s = "lambda.min", type = "response"))
pv$randomForest <- predict(m$rf, Zv[, feats, drop = FALSE], type = "prob")[, "APAP"]
pv$XGBoost      <- predict(m$xgb, xgb.DMatrix(Zv[, feats, drop = FALSE]))
write.csv(pv, file.path(out_dir, "validation_predictions.csv"), row.names = FALSE)

## ---------- 4. AUC(区分度)与校准摘要 ----------
res <- data.frame(); roc_dat <- list()
cal <- data.frame()
for (mn in c("LASSO", "randomForest", "XGBoost")) {
  rr <- roc(yv, pv[[mn]], levels = c("Control", "APAP"), direction = "<", quiet = TRUE)
  res <- rbind(res, data.frame(model = mn, auc = round(as.numeric(auc(rr)), 4)))
  cal <- rbind(cal, data.frame(
    model = mn,
    median_prob_Control = round(median(pv[[mn]][yv == "Control"]), 3),
    range_prob_Control  = sprintf("%.2f-%.2f", min(pv[[mn]][yv == "Control"]), max(pv[[mn]][yv == "Control"])),
    median_prob_APAP    = round(median(pv[[mn]][yv == "APAP"]), 3),
    range_prob_APAP     = sprintf("%.2f-%.2f", min(pv[[mn]][yv == "APAP"]), max(pv[[mn]][yv == "APAP"]))))
  co <- coords(rr, "all", ret = c("specificity", "sensitivity"), as.list = FALSE)
  roc_dat[[mn]] <- data.frame(model = mn, fpr = 1 - co$specificity,
                              tpr = co$sensitivity, auc = round(as.numeric(auc(rr)), 3))
}
write.csv(res, file.path(out_dir, "validation_auc.csv"), row.names = FALSE)
write.csv(cal, file.path(out_dir, "validation_calibration.csv"), row.names = FALSE)
print(res, row.names = FALSE); print(cal, row.names = FALSE)

roc_all <- do.call(rbind, roc_dat)
lab <- unique(roc_all[, c("model", "auc")])
p <- ggplot(roc_all, aes(fpr, tpr, color = model)) +
  geom_path(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_brewer(palette = "Set1",
                     labels = setNames(sprintf("%s (AUC=%.3f)", lab$model, lab$auc), lab$model)) +
  coord_equal() + theme_bw(base_size = 12) +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = "External validation ROC (GSE111828)",
       subtitle = "4 baseline vs 20 APAP (12-72h), z-scored with training statistics",
       color = NULL)
ggsave(file.path(out_dir, "roc_validation.png"), p, width = 6, height = 5.5, dpi = 300)
cat("[05] 外部验证完成\n")
