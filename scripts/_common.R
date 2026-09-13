# _common.R — 项目根目录、GEO 注释解析、矩阵标准化
# 由各分析脚本 source,勿单独运行

os_ali_base <- function() {
  env <- Sys.getenv("OS_ALI_BASE", unset = "")
  if (nzchar(env) && dir.exists(env)) {
    return(normalizePath(env, winslash = "/", mustWork = TRUE))
  }
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grepl("^--file=", ca)])
  if (length(f) == 1 && file.exists(f)) {
    return(normalizePath(file.path(dirname(f), ".."), winslash = "/", mustWork = TRUE))
  }
  for (cand in c(".", "..")) {
    if (file.exists(file.path(cand, "scripts")) &&
        file.exists(file.path(cand, "report"))) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }
  stop("无法确定项目根目录。请从仓库根目录运行,或设置环境变量 OS_ALI_BASE。")
}

## GEO platform annot.gz → probe, gene_symbol(取 /// 分隔的第一个符号)
parse_geo_annot <- function(path) {
  if (!file.exists(path)) stop("找不到注释文件: ", path)
  con <- gzfile(path, open = "rt")
  raw <- readLines(con, warn = FALSE)
  close(con)
  hdr <- grep("^ID[\t ]", raw)[1]
  if (is.na(hdr)) stop("注释文件无 ID 表头: ", path)
  end <- grep("^!(Annotation|platform)_table_end", raw, ignore.case = TRUE)
  last <- if (length(end)) end[1] - 1 else length(raw)
  tab <- read.delim(textConnection(raw[hdr:last]), header = TRUE,
                    sep = "\t", colClasses = "character", quote = "",
                    check.names = TRUE)
  if (!"ID" %in% names(tab)) stop("注释表缺少 ID 列: ", path)
  sym_col <- grep("gene.?symbol", names(tab), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(sym_col) || !nzchar(sym_col)) {
    stop("注释文件无 Gene symbol 列: ", paste(names(tab), collapse = ", "))
  }
  sym <- tab[[sym_col]]
  sym[is.na(sym)] <- ""
  sym <- sub(" *///.*", "", sym)
  sym <- trimws(sym)
  keep <- nzchar(sym) & !sym %in% c("---", "NA")
  out <- data.frame(probe = tab$ID[keep], gene_symbol = sym[keep],
                    stringsAsFactors = FALSE)
  out[!duplicated(out$probe), ]
}

## 按列 z-score;可传入训练集的 center/scale 应用到测试集
zscore_cols <- function(X, center = NULL, scale. = NULL) {
  X <- as.matrix(X)
  if (is.null(center)) center <- colMeans(X, na.rm = TRUE)
  if (is.null(scale.)) scale. <- apply(X, 2, stats::sd, na.rm = TRUE)
  scale.[!is.finite(scale.) | scale. == 0] <- 1
  Z <- scale(X, center = center, scale = scale.)
  Z <- matrix(as.numeric(Z), nrow = nrow(X), ncol = ncol(X), dimnames = dimnames(X))
  attr(Z, "scaled:center") <- as.numeric(center)
  names(attr(Z, "scaled:center")) <- colnames(X)
  attr(Z, "scaled:scale") <- as.numeric(scale.)
  names(attr(Z, "scaled:scale")) <- colnames(X)
  Z
}

## 分层 K 折(每折尽量保持类别比例)
stratified_folds <- function(y, k) {
  y <- as.factor(y)
  folds <- integer(length(y))
  for (lev in levels(y)) {
    ii <- which(y == lev)
    ii <- sample(ii)
    folds[ii] <- rep_len(seq_len(k), length(ii))
  }
  split(seq_along(y), folds)
}

## 将新矩阵列对齐到训练特征;缺失基因填 0(z-score 后的均值)
align_features <- function(X, feats) {
  X <- as.matrix(X)
  miss <- setdiff(feats, colnames(X))
  if (length(miss)) {
    pad <- matrix(0, nrow = nrow(X), ncol = length(miss),
                  dimnames = list(rownames(X), miss))
    X <- cbind(X, pad)
  }
  X[, feats, drop = FALSE]
}
