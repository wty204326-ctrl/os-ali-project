#!/usr/bin/env Rscript
# 04_enrichment.R — OS-DEGs 的 GO/KEGG 功能富集分析
# 用法: Rscript scripts/04_enrichment.R
# 注意:查询集本身来自 GO:0006979 交集,该条目显著是筛选的必然结果,不作独立发现。

ca <- commandArgs(trailingOnly = FALSE)
f  <- sub("^--file=", "", ca[grepl("^--file=", ca)])
script_dir <- if (length(f) == 1) dirname(normalizePath(f)) else "scripts"
source(file.path(script_dir, "_common.R"))

base_dir <- os_ali_base()
out_dir  <- file.path(base_dir, "results")
set.seed(20260909)
suppressMessages({
  library(clusterProfiler); library(org.Mm.eg.db); library(ggplot2)
})

os_deg <- read.csv(file.path(out_dir, "os_deg_GSE17649.csv"), stringsAsFactors = FALSE)
ann    <- read.csv(file.path(out_dir, "annot_GSE17649.csv"), stringsAsFactors = FALSE)

sig_genes <- unique(os_deg$gene_symbol)
universe  <- unique(ann$gene_symbol)
map1 <- suppressMessages(bitr(sig_genes, fromType = "SYMBOL",
                              toType = "ENTREZID", OrgDb = org.Mm.eg.db))
map0 <- suppressMessages(bitr(universe, fromType = "SYMBOL",
                              toType = "ENTREZID", OrgDb = org.Mm.eg.db))
cat(sprintf("OS-DEGs %d 个符号映射到 %d 个 EntrezID;背景 %d → %d\n",
            length(sig_genes), nrow(map1), length(universe), nrow(map0)))

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
cat(sprintf("GO BP 显著通路(q<0.2): %d\n", sum(res_ego$qvalue < 0.2, na.rm = TRUE)))

tryCatch({
  ekegg <- enrichKEGG(gene         = map1$ENTREZID,
                      universe     = map0$ENTREZID,
                      organism     = "mmu",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.25)
  ekegg <- setReadable(ekegg, org.Mm.eg.db, keyType = "ENTREZID")
  d <- as.data.frame(ekegg)
  write.csv(d, file.path(out_dir, "enrichment_kegg.csv"), row.names = FALSE)
  cat(sprintf("KEGG 显著通路(q<0.25): %d\n", sum(d$qvalue < 0.25, na.rm = TRUE)))
}, error = function(e) {
  cat("KEGG 在线查询失败(不影响 GO):", conditionMessage(e), "\n")
})

## 点图:去掉筛选所用的 GO:0006979,避免把循环论证画成主结果
plot_df <- res_ego[!is.na(res_ego$qvalue) & res_ego$qvalue < 0.2 & res_ego$ID != "GO:0006979", ]
top <- head(plot_df, 8)
if (nrow(top) > 0) {
  top$Description <- factor(top$Description, levels = rev(top$Description))
  p <- ggplot(top, aes(x = Count, y = Description,
                       size = Count, color = qvalue)) +
    geom_point() + scale_color_gradient(low = "#d73027", high = "#4575b4") +
    theme_bw(base_size = 12) + labs(x = "Gene count", y = NULL,
    title = "GO BP enrichment of OS-DEGs",
    subtitle = "GO:0006979 (used for gene selection) excluded")
  ggsave(file.path(out_dir, "enrichment_dotplot.png"), p,
         width = 8, height = 4.5, dpi = 300)
  cat("dotplot 已保存(已排除 GO:0006979)\n")
} else {
  cat("无显著通路可绘图\n")
}
cat("[04] 完成\n")
