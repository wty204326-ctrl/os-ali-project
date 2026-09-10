#!/usr/bin/env Rscript
# 06b_crossspecies.R — 人源 GSE120652 跨物种方向一致性检验
#
# 方法学说明(修正 v1.0 的选择性报告):
#   - 小鼠→人基因映射:7 个检验基因均为 MGI 收录的 1:1 直系同源基因,
#     符号直配(大小写转换)与 MGI HOM_MouseHumanSequence 一致;如扩展到
#     全基因组应改用 MGI 同源报告文件,此处子集规模符号直配是充分的。
#   - 人源数据 n=3v3(FFPE),无一基因过 FDR 校正,结果仅作方向性参考,
#     不支持强结论;不一致的基因如实列出。
#
# 用法: Rscript scripts/06b_crossspecies.R
# 输出: results/crossspecies_GSE120652.csv

sa <- commandArgs(FALSE)
fs <- sub("^--file=", "", sa[grepl("^--file=", sa)])
base_dir <- if (length(fs)) normalizePath(file.path(dirname(fs), "..")) else normalizePath(getwd())
out_dir <- file.path(base_dir, "results")
set.seed(20260909)
suppressMessages(library(limma))

f <- file.path(base_dir, "data/raw/GSE120652_series_matrix.txt.gz")
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
cls <- rep(NA_character_, n)
for (i in seq_len(n)) {
  vals <- paste(sapply(chars, function(cc) if (length(cc) >= i) cc[i] else NA), collapse = "|")
  if (grepl("failure|ALF|injur", vals, ignore.case = TRUE)) cls[i] <- "ALF"
  else if (grepl("normal|healthy|control", vals, ignore.case = TRUE)) cls[i] <- "Normal"
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

group <- factor(cls, levels = c("Normal", "ALF"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit2 <- eBayes(contrasts.fit(lmFit(mat, design),
                             makeContrasts(ALF - Normal, levels = design)))
de_h <- topTable(fit2, number = Inf)
de_h$probe <- rownames(de_h)

## 探针注释(GPL6244,动态表头)
ann_f <- file.path(base_dir, "data/raw/GPL6244.annot.gz")
raw <- readLines(gzfile(ann_f))
hdr <- grep("^ID\t", raw)[1]
ann_raw <- read.delim(textConnection(raw[hdr:length(raw)]), header = TRUE,
                      sep = "\t", colClasses = "character", quote = "")
ann <- ann_raw[ann_raw$ID %in% de_h$probe, c("ID", "Gene.symbol")]
names(ann) <- c("probe", "gene_symbol")
dem <- merge(de_h, ann, by = "probe")

## 小鼠 OS-DEGs → 人源直系同源(MGI HOM_MouseHumanSequence 官方映射)
os_m <- read.csv(file.path(out_dir, "os_deg_GSE17649.csv"))
mouse_genes <- unique(os_m$gene_symbol)
mgi_f <- file.path(base_dir, "data/raw/HOM_MouseHumanSequence.rpt")
if (!file.exists(mgi_f)) stop("缺少 MGI 同源报告,请先运行 01_download.sh")
mgi <- read.delim(mgi_f, check.names = FALSE, colClasses = "character")
mgi$key <- mgi[["DB Class Key"]]
ortho <- do.call(rbind, lapply(split(mgi, mgi$key), function(g) {
  mu <- g[g[["NCBI Taxon ID"]] == "10090", "Symbol"]
  hu <- g[g[["NCBI Taxon ID"]] == "9606", "Symbol"]
  if (length(mu) == 1 && length(hu) >= 1)
    data.frame(mouse = mu, human = hu, n_human = length(hu), stringsAsFactors = FALSE)
}))
ortho_m <- ortho[ortho$mouse %in% mouse_genes, ]
cat(sprintf("MGI 映射:OS-DEG %d 个小鼠基因中 %d 个有直系同源记录\n",
            length(mouse_genes), length(unique(ortho_m$mouse))))
print(ortho_m, row.names = FALSE)
human_genes <- unique(ortho_m$human)
if (length(human_genes) == 0) stop("无直系同源映射,终止")

res <- dem[dem$gene_symbol %in% human_genes, c("gene_symbol", "logFC", "P.Value", "adj.P.Val")]
res$mouse_gene <- ortho_m$mouse[match(res$gene_symbol, ortho_m$human)]
res$mouse_logFC <- sapply(res$mouse_gene, function(g) os_m$logFC[os_m$gene_symbol == g][1])
res <- res[order(-abs(res$logFC)), ]
res$direction_concordant <- sign(res$logFC) == sign(res$mouse_logFC)
write.csv(res, file.path(out_dir, "crossspecies_GSE120652.csv"), row.names = FALSE)
cat(sprintf("检验基因: %d;方向一致: %d (%s);全部 adj.P>0.15 为真时无FDR水平结论\n",
            nrow(res), sum(res$direction_concordant),
            paste(res$gene_symbol[res$direction_concordant], collapse = "/")))
print(res[, c("mouse_gene", "gene_symbol", "mouse_logFC", "logFC", "adj.P.Val",
              "direction_concordant")], row.names = FALSE)
cat("[06b] 完成(注意:结果仅供方向性参考,不支持强结论)\n")
