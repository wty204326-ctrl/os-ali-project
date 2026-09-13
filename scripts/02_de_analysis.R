#!/usr/bin/env Rscript
# 02_de_analysis.R — GEO series matrix 读取 + limma 差异表达分析
# 用法: Rscript scripts/02_de_analysis.R <GSE编号> <mouse|human>
# 输出: results/de_all_<acc>.csv      全部探针 limma 结果
#       results/deg_<acc>.csv         显著差异探针 (adj.P<0.05, |logFC|>1)
#       results/annot_<acc>.csv       由 GPL annot.gz 生成的探针注释
#       results/expr_<acc>.rds        处理后的表达矩阵与分组(供 ML 脚本复用)
#       results/os_deg_<acc>.csv      OS-DEGs

ca <- commandArgs(trailingOnly = FALSE)
f  <- sub("^--file=", "", ca[grepl("^--file=", ca)])
script_dir <- if (length(f) == 1) dirname(normalizePath(f)) else "scripts"
source(file.path(script_dir, "_common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("用法: Rscript scripts/02_de_analysis.R <GSE编号> [mouse|human]")
acc <- args[1]
species <- ifelse(length(args) >= 2, args[2], "mouse")

base_dir <- os_ali_base()
raw_dir  <- file.path(base_dir, "data/raw")
out_dir  <- file.path(base_dir, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(20260909)

suppressMessages(library(limma))

## ---------- 1. 读取 series matrix ----------
f <- file.path(raw_dir, paste0(acc, "_series_matrix.txt.gz"))
if (!file.exists(f)) stop("找不到文件: ", f)
lines <- readLines(gzfile(f), warn = FALSE)

get_meta <- function(pattern) {
  ln <- lines[grepl(pattern, lines)]
  if (length(ln) == 0) return(NULL)
  v <- unlist(strsplit(ln, "\t"))
  v <- gsub('^"|"$', "", v)
  v[v != "" & !grepl("^!", v)]
}
gsm <- get_meta("^!Sample_geo_accession")

## 分组:0hr / control / vehicle 记为 Control(GSE17649 为给药后 0h 基线,非溶媒对照)
char_lines <- lines[grepl("^!Sample_characteristics_ch1", lines)]
chars <- lapply(char_lines, function(x) {
  v <- strsplit(x, "\t")[[1]]
  gsub('^"|"$', "", v[-1])
})
n <- length(gsm)
cls <- rep(NA_character_, n)
for (i in seq_len(n)) {
  vals <- paste(sapply(chars, function(cc) if (length(cc) >= i) cc[i] else NA), collapse = "|")
  if (grepl("APAP|acetaminophen", vals, ignore.case = TRUE) &&
      !grepl("control|0hr|vehicle", vals, ignore.case = TRUE)) {
    cls[i] <- "APAP"
  } else if (grepl("control|0hr|vehicle", vals, ignore.case = TRUE)) {
    cls[i] <- "Control"
  }
}
if (any(is.na(cls))) stop("有样本无法分组: ", paste(gsm[is.na(cls)], collapse = ", "))

## ---------- 2. 表达矩阵 ----------
tb <- grep("^!series_matrix_table_begin", lines)
te <- grep("^!series_matrix_table_end", lines)
con <- textConnection(lines[(tb + 1):(te - 1)])
tab <- read.delim(con, row.names = 1, check.names = FALSE)
close(con)
mat <- as.matrix(tab)
if (!identical(colnames(mat), gsm)) {
  if (all(gsm %in% colnames(mat))) {
    mat <- mat[, gsm]
  } else {
    stop("表达矩阵列名与样本元数据不匹配")
  }
}
transformed <- FALSE
if (max(mat, na.rm = TRUE) > 100) {
  mat <- log2(mat + 1)
  transformed <- TRUE
}
cat(sprintf("[%s] samples=%d (Control=%d, APAP=%d) probes=%d log2转换=%s\n",
            acc, n, sum(cls == "Control"), sum(cls == "APAP"),
            nrow(mat), transformed))

## ---------- 3. limma 差异分析 ----------
group <- factor(cls, levels = c("Control", "APAP"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(mat, design)
cont <- makeContrasts(APAP - Control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))
de <- topTable(fit2, number = Inf)
de$probe <- rownames(de)
de <- de[, c("probe", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
write.csv(de, file.path(out_dir, paste0("de_all_", acc, ".csv")), row.names = FALSE)

deg <- subset(de, adj.P.Val < 0.05 & abs(logFC) > 1)
deg <- deg[order(deg$adj.P.Val), ]
write.csv(deg, file.path(out_dir, paste0("deg_", acc, ".csv")), row.names = FALSE)
cat(sprintf("[%s] DEGs(adj.P<0.05,|logFC|>1): %d (up=%d, down=%d)\n",
            acc, nrow(deg), sum(deg$logFC > 0), sum(deg$logFC < 0)))

## ---------- 4. 探针注释(从 series 的 GPL annot.gz 生成,不再依赖手工 CSV) ----------
plat <- get_meta("^!Series_platform_id")
plat <- if (length(plat)) plat[1] else NA_character_
ann_f <- file.path(out_dir, paste0("annot_", acc, ".csv"))
gpl_f <- if (!is.na(plat)) file.path(raw_dir, paste0(plat, ".annot.gz")) else ""
if (nzchar(gpl_f) && file.exists(gpl_f)) {
  ann <- parse_geo_annot(gpl_f)
  write.csv(ann, ann_f, row.names = FALSE)
  cat(sprintf("[%s] 注释 %s → %d 探针\n", acc, basename(gpl_f), nrow(ann)))
} else if (file.exists(ann_f)) {
  ann <- read.csv(ann_f, stringsAsFactors = FALSE)
  cat(sprintf("[%s] 使用已有注释 %s (%d 探针);未找到 %s\n",
              acc, basename(ann_f), nrow(ann), basename(gpl_f)))
} else {
  ann <- NULL
  cat(sprintf("[%s] 无探针注释(缺 %s),跳过 OS-DEGs\n", acc, basename(gpl_f)))
}

os_f <- file.path(base_dir, "data", paste0("os_genes_", species, ".csv"))
if (!is.null(ann) && file.exists(os_f)) {
  os <- read.csv(os_f, stringsAsFactors = FALSE)$gene_symbol
  dem <- merge(de, ann, by = "probe")
  os_deg <- subset(dem, gene_symbol %in% os & adj.P.Val < 0.05 & abs(logFC) > 1)
  os_deg <- os_deg[order(os_deg$adj.P.Val), ]
  write.csv(os_deg, file.path(out_dir, paste0("os_deg_", acc, ".csv")), row.names = FALSE)
  cat(sprintf("[%s] OS-DEGs: %d 探针 / %d 基因\n",
              acc, nrow(os_deg), length(unique(os_deg$gene_symbol))))
} else if (is.null(ann) || !file.exists(os_f)) {
  cat(sprintf("[%s] 跳过 OS-DEGs(缺注释或 %s)\n", acc, basename(os_f)))
}

## ---------- 5. 保存处理后的矩阵供 ML 使用 ----------
saveRDS(list(mat = mat, group = cls, acc = acc, transformed = transformed),
        file.path(out_dir, paste0("expr_", acc, ".rds")))
cat(sprintf("[%s] 完成,结果已写入 %s\n", acc, out_dir))
