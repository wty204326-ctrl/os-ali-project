#!/usr/bin/env Rscript
# 05_validation.R — GSE111828 外部数据集验证诊断模型
# 数据: RNA-seq 逐样本 counts(xlsx combined 表), 4 Control(0hr) vs 20 APAP(12-72hr)
# 用法: Rscript scripts/05_validation.R
# 输出: results/validation_predictions.csv
#       results/validation_auc.csv
#       results/roc_validation.png

base_dir <- "/home/wty204326/os-ali-project"
out_dir  <- file.path(base_dir, "results")
set.seed(20260909)
suppressMessages({
  library(readxl); library(glmnet); library(randomForest)
  library(xgboost); library(pROC); library(ggplot2)
})

## ---------- 1. 读取验证数据 ----------
f <- file.path(base_dir, "data/raw/GSE111828_results.xlsx")
d <- suppressMessages(read_excel(f, sheet = "combined"))
sample_cols <- grep("^MMR", names(d), value = TRUE)   # S01..S24
stopifnot(length(sample_cols) == 24)
sym <- d$symbol
cnt <- as.matrix(d[, sample_cols])
rownames(cnt) <- ifelse(is.na(sym), paste0("NA_", seq_along(sym)), sym)
cnt <- cnt[!duplicated(rownames(cnt)), ]

group <- c(rep("Control", 4), rep("APAP", 20))
yv <- factor(group, levels = c("Control", "APAP"))

## ---------- 2. 取 OS-DEGs 特征(log2 计数) ----------
mods <- readRDS(file.path(out_dir, "models_GSE17649.rds"))
feats <- mods$features
missing <- setdiff(feats, rownames(cnt))
cat("特征基因在验证集中缺失:", ifelse(length(missing) == 0, "无", paste(missing, collapse = ",")), "\n")
feats <- setdiff(feats, missing)
Xv <- t(log2(cnt[feats, , drop = FALSE] + 1))

## 合理性检查:对照的 Hmox1 应显著低于 APAP
if ("Hmox1" %in% feats) {
  hm <- Xv[, "Hmox1"]
  cat(sprintf("Hmox1 均值: Control=%.2f  APAP=%.2f (方向%s)\n",
              mean(hm[yv == "Control"]), mean(hm[yv == "APAP"]),
              ifelse(mean(hm[yv == "APAP"]) > mean(hm[yv == "Control"]), "正确✓", "异常✗")))
}

## ---------- 3. 三模型预测 ----------
pv <- data.frame(sample = sample_cols, group = group)
pv$LASSO        <- as.numeric(predict(mods$lasso, newx = Xv[, mods$features, drop = FALSE],
                                      s = "lambda.min", type = "response"))
pv$randomForest <- predict(mods$rf, Xv[, mods$features, drop = FALSE], type = "prob")[, "APAP"]
pv$XGBoost      <- predict(mods$xgb, xgb.DMatrix(Xv[, mods$features, drop = FALSE]))
write.csv(pv, file.path(out_dir, "validation_predictions.csv"), row.names = FALSE)

res <- data.frame()
roc_dat <- list()
for (m in c("LASSO", "randomForest", "XGBoost")) {
  rr <- roc(yv, pv[[m]], levels = c("Control", "APAP"), direction = "<", quiet = TRUE)
  res <- rbind(res, data.frame(model = m, auc = round(as.numeric(auc(rr)), 4)))
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
       subtitle = "4 Control (0hr) vs 20 APAP (12-72hr), independent RNA-seq", color = NULL)
ggsave(file.path(out_dir, "roc_validation.png"), p, width = 6, height = 5.5, dpi = 300)
cat("[05] 外部验证完成\n")
