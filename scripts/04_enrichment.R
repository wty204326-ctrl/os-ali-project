#!/usr/bin/env Rscript
# 04_enrichment.R — 全部 DEGs 的 GO/KEGG 功能富集分析
#
# 方法学说明(修正 v1.0 的循环论证):
#   v1.0 对 OS-DEGs(本身取自 GO:0006979 基因集)做富集,再报"最显著通路
#   为 response to oxidative stress"属循环论证。本版对**全部差异表达基因**
#   做富集,得到独立于基因集定义的生物学图景;OS-DEGs 与 GO:0006979 的
#   交集是定义使然,不再作为富集证据呈现。
#   显著性阈值:GO/KEGG 均为 p.adjust < 0.05(正文只引用 q<0.05 的条目)。
#
# 用法: Rscript scripts/04_enrichment.R
# 输出: results/enrichment_go_BP.csv
#       results/enrichment_kegg.csv(网络可达时)

sa <- commandArgs(FALSE)
fs <- sub("^--file=", "", sa[grepl("^--file=", sa)])
base_dir <- if (length(fs)) normalizePath(file.path(dirname(fs), "..")) else normalizePath(getwd())
out_dir <- file.path(base_dir, "results")
set.seed(20260909)
suppressMessages({
  library(clusterProfiler); library(org.Mm.eg.db); library(ggplot2)
})

deg <- read.csv(file.path(out_dir, "deg_GSE17649.csv"))
ann <- read.csv(file.path(out_dir, "annot_GSE17649.csv"))
sig_genes <- unique(ann$gene_symbol[ann$probe %in% deg$probe])
universe <- unique(ann$gene_symbol)
m1 <- suppressMessages(bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
m0 <- suppressMessages(bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
cat(sprintf("全部 DEGs %d 符号映射 %d EntrezID;背景 %d → %d\n",
            length(sig_genes), nrow(m1), length(universe), nrow(m0)))

## GO BP(全 DEGs)
ego <- enrichGO(gene = m1$ENTREZID, universe = m0$ENTREZID,
                OrgDb = org.Mm.eg.db, ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1,
                readable = TRUE)
res_ego <- as.data.frame(ego)
write.csv(res_ego, file.path(out_dir, "enrichment_go_BP.csv"), row.names = FALSE)
cat(sprintf("GO BP 显著通路(p.adjust<0.05): %d\n", sum(res_ego$p.adjust < 0.05)))

## KEGG(全 DEGs;在线查询,失败不阻塞)
res_k <- tryCatch({
  ek <- enrichKEGG(gene = m1$ENTREZID, universe = m0$ENTREZID, organism = "mmu",
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1)
  ek <- setReadable(ek, org.Mm.eg.db, keyType = "ENTREZID")
  d <- as.data.frame(ek)
  write.csv(d, file.path(out_dir, "enrichment_kegg.csv"), row.names = FALSE)
  cat(sprintf("KEGG 显著通路(p.adjust<0.05): %d\n", sum(d$p.adjust < 0.05)))
  d
}, error = function(e) {
  cat("KEGG 在线查询失败(不影响 GO):", conditionMessage(e), "\n"); NULL
})

## 点图(GO 前 8)
top <- head(res_ego[order(res_ego$p.adjust), ], 8)
if (nrow(top) > 0) {
  top$Description <- factor(top$Description, levels = rev(top$Description))
  p <- ggplot(top, aes(x = Count, y = Description, size = Count, color = p.adjust)) +
    geom_point() + scale_color_gradient(low = "#d73027", high = "#4575b4") +
    theme_bw(base_size = 12) +
    labs(x = "gene count", y = NULL, title = "GO BP enrichment of all DEGs (top 8)")
  ggsave(file.path(out_dir, "enrichment_dotplot.png"), p, width = 8.5, height = 4.5, dpi = 300)
  cat("dotplot 已保存\n")
}
cat("[04] 完成\n")
