#!/bin/bash
# run_all.sh — 一键复现全部分析
# 用法: bash scripts/run_all.sh
set -e
BASE="$(cd "$(dirname "$0")/.." && pwd)"
source ~/miniforge3/bin/activate rna

echo "===== [0/6] 依赖自检与补装 ====="
bash "$BASE/scripts/00_install_deps.sh"

echo "===== [1/6] 数据下载 ====="
bash "$BASE/scripts/01_download.sh"

cd "$BASE"
echo "===== [2/6] 差异表达分析(训练集 GSE17649)====="
Rscript scripts/02_de_analysis.R GSE17649 mouse

echo "===== [3/6] 机器学习建模 ====="
Rscript scripts/03_ml_modeling.R

echo "===== [4/6] 功能富集 ====="
Rscript scripts/04_enrichment.R

echo "===== [5/6] 外部验证(GSE111828)====="
Rscript scripts/05_validation.R

echo "===== [6/6] 跨物种方向验证 + 报告渲染 ====="
Rscript scripts/06b_crossspecies.R
Rscript -e "rmarkdown::render('report/06_report.Rmd', quiet = TRUE)"

echo ""
echo "全部完成 ✔  课程报告: report/06_report.html"
