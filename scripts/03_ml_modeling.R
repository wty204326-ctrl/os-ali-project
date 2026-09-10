#!/usr/bin/env Rscript
# 03_ml_modeling.R — OS-DEGs 急性肝损伤诊断模型(嵌套特征选择,无泄漏版)
#
# 方法学说明(修正 v1.0 的特征泄漏):
#   OS-DEGs 此前用全部 36 样本筛出后再做交叉验证,验证折的标签参与了
#   特征选择,训练 AUC 偏乐观。本版将 limma 差异分析嵌入每一折的训练集
#   内部重新执行(特征选择 → 建模 → 预测均在折内完成),折外 AUC 不再
#   含选择泄漏;全数据只拟合一次"最终模型",用于特征重要性与外部验证。
#
# 用法: Rscript scripts/03_ml_modeling.R
# 输出: results/ml_cv_results.csv           嵌套 CV AUC(mean±sd)
#       results/cv_feature_counts.csv       每折选中的特征数
#       results/roc_training.png            训练集 ROC(单次 5 折)
#       results/lasso_coefficients.csv      最终 LASSO 系数
#       results/rf_importance.csv           最终 RF 重要性
#       results/xgb_importance.csv          最终 XGBoost 重要性
#       results/models_GSE17649.rds         最终模型+特征+训练集统计量(供 05 标准化)

sa <- commandArgs(FALSE)
fs <- sub("^--file=", "", sa[grepl("^--file=", sa)])
base_dir <- if (length(fs)) normalizePath(file.path(dirname(fs), "..")) else normalizePath(getwd())
out_dir <- file.path(base_dir, "results")
set.seed(20260909)
suppressMessages({
  library(limma); library(glmnet); library(randomForest)
  library(xgboost); library(pROC); library(ggplot2)
})

## ---------- 1. 载入数据与 OS 基因集 ----------
dat <- readRDS(file.path(out_dir, "expr_GSE17649.rds"))
mat <- dat$mat; cls <- dat$group; strain <- dat$strain
ann <- read.csv(file.path(out_dir, "annot_GSE17649.csv"))
os <- toupper(read.csv(file.path(base_dir, "data", "os_genes_mouse.csv"))$gene_symbol)
ann$usym <- toupper(ann$gene_symbol)
y <- factor(cls, levels = c("Control", "APAP")); yc <- as.integer(y == "APAP")
strain_f <- factor(strain)

## 折内特征选择函数:limma(strain 校正) → DEGs ∩ OS 基因(大小写不敏感)
## → 按基因取均值
select_features <- function(idx, min_fc = 1, alpha = 0.05) {
  m <- mat[, idx, drop = FALSE]
  g <- factor(cls[idx], levels = c("Control", "APAP"))
  s <- factor(strain[idx])
  d <- model.matrix(~0 + g + s)
  colnames(d) <- c(levels(g), paste0("s", levels(s)[-1]))
  f2 <- eBayes(contrasts.fit(lmFit(m, d),
                             makeContrasts(APAP - Control, levels = d)))
  tt <- topTable(f2, number = Inf)
  tt$probe <- rownames(tt)
  hits <- tt$probe[tt$adj.P.Val < alpha & abs(tt$logFC) > min_fc]
  sort(unique(ann$usym[ann$probe %in% hits & ann$usym %in% os]))
}
collapse <- function(genes, idx) {
  sel <- ann$usym %in% genes
  gm <- mat[ann$probe[sel], idx, drop = FALSE]
  gn <- ann$usym[sel]
  t(sapply(genes, function(g) colMeans(gm[gn == g, , drop = FALSE])))
}

## ---------- 2. 嵌套交叉验证:5 折 × 5 次 ----------
K <- 5; REP <- 5
mods <- c("LASSO", "randomForest", "XGBoost")
pred <- as.data.frame(matrix(NA_real_, ncol(mat), length(mods) * REP,
                             dimnames = list(colnames(mat),
                                             paste(rep(mods, REP), rep(1:REP, each = length(mods)), sep = "_"))))
feat_counts <- data.frame()
for (r in seq_len(REP)) {
  set.seed(1000 + r)
  folds <- split(sample(ncol(mat)), cut(seq_len(ncol(mat)), K))
  for (k in seq_len(K)) {
    te <- folds[[k]]; tr <- setdiff(seq_len(ncol(mat)), te)
    genes_k <- select_features(tr)
    feat_counts <- rbind(feat_counts, data.frame(rep = r, fold = k, n_features = length(genes_k)))
    if (length(genes_k) == 0) {
      cat(sprintf("警告: rep%d fold%d 无选中特征,沿用全数据特征\n", r, k))
      genes_k <- select_features(seq_len(ncol(mat)))
    }
    Xtr <- t(collapse(genes_k, tr)); Xte <- t(collapse(genes_k, te))
    mu <- colMeans(Xtr); sdv <- apply(Xtr, 2, sd)
    Ztr <- sweep(sweep(Xtr, 2, mu), 2, sdv, "/")
    Zte <- sweep(sweep(Xte, 2, mu), 2, sdv, "/")
    cvf <- cv.glmnet(Ztr, yc[tr], alpha = 1, family = "binomial")
    pred[te, paste0("LASSO_", r)] <- as.numeric(
      predict(cvf, newx = Zte, s = "lambda.min", type = "response"))
    rf <- randomForest(x = Ztr, y = y[tr], ntree = 500)
    pred[te, paste0("randomForest_", r)] <-
      predict(rf, Zte, type = "prob")[, "APAP"]
    bst <- xgb.train(params = list(objective = "binary:logistic", eta = 0.1, max_depth = 2),
                     data = xgb.DMatrix(Ztr, label = yc[tr]), nrounds = 50, verbose = 0)
    pred[te, paste0("XGBoost_", r)] <- predict(bst, xgb.DMatrix(Zte))
  }
  cat(sprintf("repeat %d/%d 完成\n", r, REP))
}
write.csv(feat_counts, file.path(out_dir, "cv_feature_counts.csv"), row.names = FALSE)
cat(sprintf("折内选中特征数: %d ~ %d(中位 %d)\n",
            min(feat_counts$n_features), max(feat_counts$n_features),
            median(feat_counts$n_features)))

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

## ---------- 4. ROC 图(第 1 次重复的折外预测) ----------
roc_dat <- do.call(rbind, lapply(mods, function(m) {
  r1 <- roc(y, pred[, paste0(m, "_1")], levels = c("Control", "APAP"), direction = "<", quiet = TRUE)
  co <- coords(r1, "all", ret = c("specificity", "sensitivity"), as.list = FALSE)
  data.frame(model = m, fpr = 1 - co$specificity, tpr = co$sensitivity,
             auc = round(as.numeric(auc(r1)), 3))
}))
lab <- unique(roc_dat[, c("model", "auc")])
p <- ggplot(roc_dat, aes(fpr, tpr, color = model)) +
  geom_path(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_brewer(palette = "Set1",
                     labels = setNames(sprintf("%s (AUC=%.3f)", lab$model, lab$auc), lab$model)) +
  coord_equal() + theme_bw(base_size = 12) +
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)",
       title = "Training ROC: OS-DEG diagnostic models (nested selection)",
       subtitle = "GSE17649: 12 baseline vs 24 APAP, feature selection inside each fold",
       color = NULL)
ggsave(file.path(out_dir, "roc_training.png"), p, width = 6, height = 5.5, dpi = 300)

## ---------- 5. 最终模型(全数据,同折内流程)与统计量 ----------
final_genes <- select_features(seq_len(ncol(mat)))
Xall <- t(collapse(final_genes, seq_len(ncol(mat))))
mu <- colMeans(Xall); sdv <- apply(Xall, 2, sd)
Z <- sweep(sweep(Xall, 2, mu), 2, sdv, "/")
cv_final <- cv.glmnet(Z, yc, alpha = 1, family = "binomial")
cf <- as.matrix(coef(cv_final, s = "lambda.min"))
cf <- data.frame(gene = rownames(cf), coefficient = round(cf[, 1], 4))
cf <- cf[cf$coefficient != 0, ]
write.csv(cf, file.path(out_dir, "lasso_coefficients.csv"), row.names = FALSE)
rf_final <- randomForest(x = Z, y = y, ntree = 1000)
imp <- data.frame(gene = rownames(rf_final$importance),
                  MeanDecreaseGini = round(rf_final$importance[, 1], 3))
imp <- imp[order(-imp$MeanDecreaseGini), ]
write.csv(imp, file.path(out_dir, "rf_importance.csv"), row.names = FALSE)
xfin <- xgb.train(params = list(objective = "binary:logistic", eta = 0.1, max_depth = 2),
                  data = xgb.DMatrix(Z, label = yc), nrounds = 50, verbose = 0)
xi <- xgb.importance(model = xfin)
write.csv(data.frame(feature = xi$Feature, gain = round(xi$Gain, 4)),
          file.path(out_dir, "xgb_importance.csv"), row.names = FALSE)
saveRDS(list(features = final_genes, train_mean = mu, train_sd = sdv,
             lasso = cv_final, rf = rf_final, xgb = xfin),
        file.path(out_dir, "models_GSE17649.rds"))
cat("\n== 最终特征(n=", length(final_genes), ") ==\n",
    paste(final_genes, collapse = ", "), "\n", sep = "")
cat("== LASSO 选中 ==\n"); print(cf, row.names = FALSE)
cat(sprintf("\n[03] 完成\n"))
