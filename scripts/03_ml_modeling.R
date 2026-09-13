#!/usr/bin/env Rscript
# 03_ml_modeling.R — 基于 OS-DEGs 的急性肝损伤机器学习诊断模型
# 模型: LASSO / randomForest / XGBoost
# 评估: 分层 5 折 × 5 次重复(折外 AUC);特征在训练集内 z-score,测试折用训练集均值/标准差
# 用法: Rscript scripts/03_ml_modeling.R

ca <- commandArgs(trailingOnly = FALSE)
f  <- sub("^--file=", "", ca[grepl("^--file=", ca)])
script_dir <- if (length(f) == 1) dirname(normalizePath(f)) else "scripts"
source(file.path(script_dir, "_common.R"))

base_dir <- os_ali_base()
out_dir  <- file.path(base_dir, "results")
set.seed(20260909)

suppressMessages({
  library(glmnet); library(randomForest); library(xgboost)
  library(pROC); library(ggplot2)
})

## ---------- 1. 特征矩阵:OS-DEGs 多探针按基因取均值 ----------
d   <- read.csv(file.path(out_dir, "os_deg_GSE17649.csv"), stringsAsFactors = FALSE)
dat <- readRDS(file.path(out_dir, "expr_GSE17649.rds"))
mat <- dat$mat; cls <- dat$group
gm  <- mat[d$probe, , drop = FALSE]
gene_uni <- unique(d$gene_symbol)
feat_mat <- t(sapply(gene_uni, function(g) colMeans(gm[d$gene_symbol == g, , drop = FALSE])))
X <- t(feat_mat)
y  <- factor(cls, levels = c("Control", "APAP"))
yc <- as.integer(y == "APAP")
cat(sprintf("特征=%d 个OS-DEGs; 样本=%d (Control=%d, APAP=%d)\n",
            ncol(X), nrow(X), sum(yc == 0), sum(yc == 1)))
cat("说明: OS-DEGs 在全训练集上选出,CV AUC 略偏乐观;以外部验证为确认指标。\n")

## ---------- 2. 分层交叉验证:5 折 × 5 次;折内 z-score ----------
K <- 5; REP <- 5
mods <- c("LASSO", "randomForest", "XGBoost")
pred <- as.data.frame(matrix(NA_real_, nrow(X), length(mods) * REP,
                             dimnames = list(rownames(X),
                                             paste(rep(mods, REP), rep(1:REP, each = length(mods)), sep = "_"))))
for (r in seq_len(REP)) {
  set.seed(1000 + r)
  folds <- stratified_folds(y, K)
  for (k in seq_along(folds)) {
    te <- folds[[k]]
    tr <- setdiff(seq_len(nrow(X)), te)
    Xtr <- zscore_cols(X[tr, , drop = FALSE])
    Xte <- zscore_cols(X[te, , drop = FALSE],
                       center = attr(Xtr, "scaled:center"),
                       scale. = attr(Xtr, "scaled:scale"))
    ## LASSO
    cvfit <- cv.glmnet(Xtr, yc[tr], alpha = 1, family = "binomial")
    pred[te, paste0("LASSO_", r)] <- as.numeric(
      predict(cvfit, newx = Xte, s = "lambda.min", type = "response"))
    ## randomForest
    rf <- randomForest(x = Xtr, y = y[tr], ntree = 500)
    pred[te, paste0("randomForest_", r)] <-
      predict(rf, Xte, type = "prob")[, "APAP"]
    ## XGBoost
    dtr <- xgb.DMatrix(Xtr, label = yc[tr])
    bst <- xgb.train(params = list(objective = "binary:logistic",
                                   eta = 0.1, max_depth = 2),
                     data = dtr, nrounds = 50, verbose = 0)
    pred[te, paste0("XGBoost_", r)] <- predict(bst, xgb.DMatrix(Xte))
  }
  cat(sprintf("repeat %d/%d 完成\n", r, REP))
}

## ---------- 3. AUC 汇总 ----------
res <- data.frame()
for (m in mods) {
  aucs <- sapply(seq_len(REP), function(r) {
    as.numeric(auc(roc(y, pred[, paste0(m, "_", r)],
                       levels = c("Control", "APAP"), direction = "<"), quiet = TRUE))
  })
  res <- rbind(res, data.frame(model = m,
                               auc_mean = round(mean(aucs), 4),
                               auc_sd = round(sd(aucs), 4)))
}
write.csv(res, file.path(out_dir, "ml_cv_results.csv"), row.names = FALSE)
print(res, row.names = FALSE)

## ---------- 4. ROC 曲线图(第 1 次重复的折外预测) ----------
roc_dat <- do.call(rbind, lapply(mods, function(m) {
  r1 <- roc(y, pred[, paste0(m, "_1")], levels = c("Control", "APAP"), direction = "<")
  co <- coords(r1, "all", ret = c("specificity", "sensitivity"), as.list = FALSE)
  data.frame(model = m, fpr = 1 - co$specificity, tpr = co$sensitivity,
             auc = round(as.numeric(auc(r1)), 3))
}))
roc_lab <- unique(roc_dat[, c("model", "auc")])
roc_lab$label <- sprintf("%s (AUC=%.3f)", roc_lab$model, roc_lab$auc)
p <- ggplot(roc_dat, aes(fpr, tpr, color = model)) +
  geom_path(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_brewer(palette = "Set1", labels = setNames(roc_lab$label, roc_lab$model)) +
  coord_equal() + theme_bw(base_size = 12) +
  labs(x = "1 − Specificity (FPR)", y = "Sensitivity (TPR)",
       title = "Training ROC: OS-DEG diagnostic models",
       subtitle = "GSE17649: 12 Control vs 24 APAP, stratified 5-fold CV (z-scored)",
       color = NULL)
ggsave(file.path(out_dir, "roc_training.png"), p, width = 6, height = 5.5, dpi = 300)

## ---------- 5. 最终模型(全数据 z-score 后拟合) ----------
Xz <- zscore_cols(X)
cv_final <- cv.glmnet(Xz, yc, alpha = 1, family = "binomial")
cf <- as.matrix(coef(cv_final, s = "lambda.min"))
cf <- data.frame(gene = rownames(cf), coefficient = round(cf[, 1], 4))
cf <- cf[cf$coefficient != 0, , drop = FALSE]
write.csv(cf, file.path(out_dir, "lasso_coefficients.csv"), row.names = FALSE)

rf_final <- randomForest(x = Xz, y = y, ntree = 1000)
imp <- data.frame(gene = rownames(rf_final$importance),
                  MeanDecreaseGini = round(rf_final$importance[, 1], 3))
imp <- imp[order(-imp$MeanDecreaseGini), ]
write.csv(imp, file.path(out_dir, "rf_importance.csv"), row.names = FALSE)

dall <- xgb.DMatrix(Xz, label = yc)
xfin <- xgb.train(params = list(objective = "binary:logistic", eta = 0.1, max_depth = 2),
                  data = dall, nrounds = 50, verbose = 0)
xi <- xgb.importance(model = xfin)
xi <- data.frame(feature = xi$Feature, gain = round(xi$Gain, 4))
write.csv(xi, file.path(out_dir, "xgb_importance.csv"), row.names = FALSE)

saveRDS(list(features = gene_uni,
             lasso = cv_final, rf = rf_final, xgb = xfin,
             scaled = TRUE),
        file.path(out_dir, "models_GSE17649.rds"))

cat("\n== LASSO 选中基因 ==\n"); print(cf, row.names = FALSE)
cat("\n== RF 重要性 Top3 ==\n"); print(head(imp, 3), row.names = FALSE)
cat(sprintf("\n[03] 完成,图与表已写入 %s\n", out_dir))
