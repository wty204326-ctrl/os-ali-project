#!/usr/bin/env Rscript
# 06b_crossspecies.R — 人源 GSE120652 跨物种方向验证
# 将小鼠 OS-DEGs 映射到人源同源基因(符号直配),在 APAP-ALF vs 正常肝 中检验 logFC 方向一致性
# 用法: Rscript scripts/06b_crossspecies.R
# 输出: results/crossspecies_GSE120652.csv

base_dir <- "/home/wty204326/os-ali-project"
out_dir  <- file.path(base_dir, "results")
set.seed(20260909)
suppressMessages(library(limma))

## ---------- 1. 读取人源 series matrix ----------
f <- file.path(base_dir, "data/raw/GSE120652_series_matrix.txt.gz")
lines <- readLines(gzfile(f))
get_meta <- function(pattern) {
  ln <- lines[grepl(pattern, lines)]
  if (length(ln) == 0) return(NULL)
  v <- unlist(strsplit(ln, "\t"))
  v <- gsub('^"|"$', "", v)
  v[v != "" & !grepl("^!", v)]
}
gsm <- get_meta("^!Sample_geo_accession")
char_lines <- lines[grepl("^!Sample_characteristics_ch1", lines)]
chars <- lapply(char_lines, function(x) {
  v <- gsub('^"|"$', "", strsplit(x, "\t")[[1]])
  v[-1]
})
n <- length(gsm)
cls <- rep(NA_character_, n)
for (i in seq_len(n)) {
  vals <- paste(sapply(chars, function(cc) if (length(cc) >= i) cc[i] else NA), collapse = "|")
  if (grepl("failure|ALF|injur", vals, ignore.case = TRUE)) {
    cls[i] <- "ALF"
  } else if (grepl("normal|healthy|control", vals, ignore.case = TRUE)) {
    cls[i] <- "Normal"
  }
}
stopifnot(!any(is.na(cls)))

tb <- grep("^!series_matrix_table_begin", lines)
te <- grep("^!series_matrix_table_end", lines)
con <- textConnection(lines[(tb + 1):(te - 1)])
tab <- read.delim(con, row.names = 1, check.names = FALSE)
close(con)
mat <- as.matrix(tab)
if (!identical(colnames(mat), gsm)) mat <- mat[, gsm]
if (max(mat, na.rm = TRUE) > 100) mat <- log2(mat + 1)
cat(sprintf("[GSE120652] samples=%d (Normal=%d, ALF=%d) probes=%d\n",
            n, sum(cls == "Normal"), sum(cls == "ALF"), nrow(mat)))

## ---------- 2. limma:ALF vs Normal ----------
group <- factor(cls, levels = c("Normal", "ALF"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit2 <- eBayes(contrasts.fit(lmFit(mat, design),
                             makeContrasts(ALF - Normal, levels = design)))
de_h <- topTable(fit2, number = Inf)
de_h$probe <- rownames(de_h)

## ---------- 3. 探针注释(GPL6244 annot,列"Gene symbol") ----------
ann_f <- file.path(base_dir, "data/raw/GPL6244.annot.gz")
if (!file.exists(ann_f)) stop("缺少 GPL6244 注释文件,请先下载")
raw <- readLines(gzfile(ann_f))
hdr <- grep("^ID\t", raw)[1]
ann_raw <- read.delim(textConnection(raw[hdr:length(raw)]), header = TRUE,
                      sep = "\t", colClasses = "character", quote = "")
ann <- ann_raw[ann_raw$ID %in% de_h$probe, c("ID", "Gene.symbol")]
names(ann) <- c("probe", "gene_symbol")
ann <- ann[!is.na(ann$gene_symbol) & ann$gene_symbol != "", ]
dem <- merge(de_h, ann, by = "probe")

## ---------- 4. 与小鼠 OS-DEGs 的方向一致性 ----------
os_m <- read.csv(file.path(out_dir, "os_deg_GSE17649.csv"))
mouse_genes <- unique(toupper(os_m$gene_symbol))          # 小鼠符号大写即人源直配
res <- dem[dem$gene_symbol %in% mouse_genes, c("gene_symbol", "logFC", "P.Value", "adj.P.Val")]
res$mouse_logFC <- sapply(res$gene_symbol, function(g) {
  os_m$logFC[toupper(os_m$gene_symbol) == g][1]
})
res <- res[order(-abs(res$logFC)), ]
res$direction_concordant <- sign(res$logFC) == sign(res$mouse_logFC)
write.csv(res, file.path(out_dir, "crossspecies_GSE120652.csv"), row.names = FALSE)
cat(sprintf("跨物种检验基因数: %d;方向一致: %d\n",
            nrow(res), sum(res$direction_concordant)))
print(res, row.names = FALSE)
cat("[06b] 完成\n")
