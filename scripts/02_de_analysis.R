#!/usr/bin/env Rscript
# 02_de_analysis.R — GSE17649 预处理 + strain 校正的 limma 差异表达分析
# 用法: Rscript scripts/02_de_analysis.R <GSE编号> <mouse|human>
# 设计说明:
#   - 分组:APAP(3h/6h 处理) vs 0h 基线对照;同一实验的 0h 并非独立健康肝,
#     该对照定义的局限在报告"局限性"中说明
#   - 品系(SJL/SMJ/DBA/C57)作为协变量纳入设计矩阵,消除品系混杂
#   - 探针注释(GPL1261→基因符号)在本脚本内解析,不再依赖历史结果文件
# 输出: results/de_all_<acc>.csv        strain 校正后的 limma 全表
#       results/deg_<acc>.csv           显著差异探针
#       results/annot_<acc>.csv         探针-基因符号映射
#       results/os_deg_<acc>.csv        OS-DEGs(与 GO:0006979 基因集交集)
#       results/expr_<acc>.rds          表达矩阵+分组+品系(供下游脚本)

sa <- commandArgs(FALSE)
fs <- sub("^--file=", "", sa[grepl("^--file=", sa)])
base_dir <- if (length(fs)) normalizePath(file.path(dirname(fs), "..")) else normalizePath(getwd())
targs <- commandArgs(trailingOnly = TRUE)
acc <- targs[1]
species <- ifelse(length(targs) >= 2, targs[2], "mouse")

raw_dir <- file.path(base_dir, "data", "raw")
out_dir <- file.path(base_dir, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(20260909)
suppressMessages(library(limma))

## ---------- 1. 读取 series matrix ----------
f <- file.path(raw_dir, paste0(acc, "_series_matrix.txt.gz"))
if (!file.exists(f)) stop("找不到文件: ", f)
lines <- readLines(gzfile(f))
get_meta <- function(pattern) {
  ln <- lines[grepl(pattern, lines)]
  if (length(ln) == 0) return(NULL)
  v <- gsub('^"|"$', "", unlist(strsplit(ln, "\t")))
  v[v != "" & !grepl("^!", v)]
}
gsm <- get_meta("^!Sample_geo_accession")
char_lines <- lines[grepl("^!Sample_characteristics_ch1", lines)]
chars <- lapply(char_lines, function(x) {
  v <- gsub('^"|"$', "", strsplit(x, "\t")[[1]]); v[-1]
})
n <- length(gsm)
cls <- rep(NA_character_, n); strain <- rep(NA_character_, n)
for (i in seq_len(n)) {
  vals <- sapply(chars, function(cc) if (length(cc) >= i) cc[i] else NA)
  joined <- paste(vals, collapse = "|")
  if (grepl("APAP|acetaminophen", joined, ignore.case = TRUE) &&
      !grepl("control|0hr|vehicle", joined, ignore.case = TRUE)) {
    cls[i] <- "APAP"
  } else if (grepl("control|0hr|vehicle", joined, ignore.case = TRUE)) {
    cls[i] <- "Control"
  }
  s <- grep("strain", vals, ignore.case = TRUE, value = TRUE)
  if (length(s)) strain[i] <- trimws(sub("^strain:", "", s[1]))
}
stopifnot(!any(is.na(cls)))

## ---------- 2. 表达矩阵 ----------
tb <- grep("^!series_matrix_table_begin", lines)
te <- grep("^!series_matrix_table_end", lines)
con <- textConnection(lines[(tb + 1):(te - 1)])
tab <- read.delim(con, row.names = 1, check.names = FALSE)
close(con)
mat <- as.matrix(tab)
if (!identical(colnames(mat), gsm)) {
  if (all(gsm %in% colnames(mat))) mat <- mat[, gsm] else stop("表达矩阵列名与样本元数据不匹配")
}
transformed <- FALSE
if (max(mat, na.rm = TRUE) > 100) { mat <- log2(mat + 1); transformed <- TRUE }
cat(sprintf("[%s] samples=%d (Control=%d, APAP=%d) probes=%d log2=%s 品系=%s\n",
            acc, n, sum(cls == "Control"), sum(cls == "APAP"), nrow(mat),
            transformed, paste(unique(na.omit(strain)), collapse = "/")))

## ---------- 3. strain 校正的 limma ----------
group <- factor(cls, levels = c("Control", "APAP"))
strain_f <- factor(strain)
design <- model.matrix(~0 + group + strain_f)
colnames(design) <- c(levels(group), paste0("strain_", levels(strain_f)[-1]))
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
cat(sprintf("[%s] DEGs(strain校正): %d (up=%d, down=%d)\n",
            acc, nrow(deg), sum(deg$logFC > 0), sum(deg$logFC < 0)))

## ---------- 4. 探针注释(在本脚本内解析 GPL1261.annot.gz) ----------
ann_f <- file.path(raw_dir, "GPL1261.annot.gz")
if (!file.exists(ann_f)) stop("缺少注释文件,请先运行 01_download.sh")
raw <- readLines(gzfile(ann_f))
hdr <- grep("^ID\t", raw)[1]
ann_raw <- read.delim(textConnection(raw[hdr:length(raw)]), header = TRUE,
                      sep = "\t", colClasses = "character", quote = "")
ann <- ann_raw[ann_raw$ID %in% rownames(mat), c("ID", "Gene.symbol")]
names(ann) <- c("probe", "gene_symbol")
ann <- subset(ann, !is.na(gene_symbol) & gene_symbol != "")
write.csv(ann, file.path(out_dir, paste0("annot_", acc, ".csv")), row.names = FALSE)
cat(sprintf("[%s] 有基因符号注释的探针: %d/%d (%.0f%%);无符号探针无法参与基因集交集\n",
            acc, nrow(ann), nrow(mat), 100 * nrow(ann) / nrow(mat)))

## ---------- 5. OS-DEGs(符号大小写不敏感匹配) ----------
os_f <- file.path(base_dir, "data", paste0("os_genes_", species, ".csv"))
if (file.exists(os_f)) {
  os <- toupper(read.csv(os_f)$gene_symbol)
  dem <- merge(de, ann, by = "probe")
  dem$up_sym <- toupper(dem$gene_symbol)
  os_deg <- subset(dem, up_sym %in% os & adj.P.Val < 0.05 & abs(logFC) > 1)
  os_deg <- os_deg[order(os_deg$adj.P.Val), ]
  os_deg$up_sym <- NULL
  write.csv(os_deg, file.path(out_dir, paste0("os_deg_", acc, ".csv")), row.names = FALSE)
  cat(sprintf("[%s] OS-DEGs: %d 探针 / %d 基因\n",
              acc, nrow(os_deg), length(unique(os_deg$gene_symbol))))
} else {
  cat(sprintf("[%s] 跳过 OS-DEGs(缺 %s)\n", acc, basename(os_f)))
}

## ---------- 6. 保存供下游使用 ----------
saveRDS(list(mat = mat, group = cls, strain = strain, acc = acc,
             transformed = transformed),
        file.path(out_dir, paste0("expr_", acc, ".rds")))
cat(sprintf("[%s] 完成\n", acc))
