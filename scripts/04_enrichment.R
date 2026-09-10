#!/usr/bin/env Rscript
# 04_enrichment.R — OS-DEGs 的 GO/KEGG 功能富集分析
# 用法: Rscript scripts/04_enrichment.R
# 输出: results/enrichment_go_BP.csv
#       results/enrichment_kegg.csv(网络可达时)
#       results/enrichment_dotplot.png

base_dir <- "/home/wty204326/os-ali-project"
out_dir  <- file.path(base_dir, "results")
set.seed(20260909)
suppressMessages({
  library(clusterProfiler); library(org.Mm.eg.db); library(ggplot2)
})

os_deg <- read.csv(file.path(out_dir, "os_deg_GSE17649.csv"))
ann    <- read.csv(file.path(out_dir, "annot_GSE17649.csv"))

## 差异基因与背景基因集:符号 → Entrez ID
sig_genes <- unique(os_deg$gene_symbol)
universe  <- unique(ann$gene_symbol)
map1 <- suppressMessages(bitr(sig_genes, fromType = "SYMBOL",
                              toType = "ENTREZID", OrgDb = org.Mm.eg.db))
map0 <- suppressMessages(bitr(universe, fromType = "SYMBOL",
                              toType = "ENTREZID", OrgDb = org.Mm.eg.db))
cat(sprintf("OS-DEGs %d 个符号映射到 %d 个 EntrezID;背景 %d → %d\n",
            length(sig_genes), nrow(map1), length(universe), nrow(map0)))

## GO BP 富集
ego <- enrichGO(gene          = map1$ENTREZID,
                universe      = map0$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.2,
                readable      = TRUE)
res_ego <- as.data.frame(ego)
write.csv(res_ego, file.path(out_dir, "enrichment_go_BP.csv"), row.names = FALSE)
cat(sprintf("GO BP 显著通路(q<0.2): %d\n", sum(res_ego$qvalue < 0.2)))

## KEGG 富集(在线查询,失败不阻塞)
res_k <- tryCatch({
  ekegg <- enrichKEGG(gene         = map1$ENTREZID,
                      universe     = map0$ENTREZID,
                      organism     = "mmu",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.25)
  ekegg <- setReadable(ekegg, org.Mm.eg.db, keyType = "ENTREZID")
  d <- as.data.frame(ekegg)
  write.csv(d, file.path(out_dir, "enrichment_kegg.csv"), row.names = FALSE)
  cat(sprintf("KEGG 显著通路(q<0.25): %d\n", sum(d$qvalue < 0.25)))
  d
}, error = function(e) {
  cat("KEGG 在线查询失败(不影响 GO):", conditionMessage(e), "\n")
  NULL
})

## 点图
top <- rbind(head(res_ego[res_ego$qvalue < 0.2, ], 8))
if (nrow(top) > 0) {
  top$Description <- factor(top$Description, levels = rev(top$Description))
  p <- ggplot(top, aes(x = Count, y = Description,
                       size = Count, color = qvalue)) +
    geom_point() + scale_color_gradient(low = "#d73027", high = "#4575b4") +
    theme_bw(base_size = 12) + labs(x = "基因数", y = NULL,
    title = "OS-DEGs 的 GO 生物过程富集(前8)")
  ggsave(file.path(out_dir, "enrichment_dotplot.png"), p,
         width = 8, height = 4.5, dpi = 300)
  cat("dotplot 已保存\n")
} else {
  cat("无显著通路可绘图\n")
}
cat("[04] 完成\n")
