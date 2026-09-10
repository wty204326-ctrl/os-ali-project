#!/bin/bash
# 01_download.sh — 数据与基因集下载(GEO + EBI GOA)
# 每个文件带完整性校验(gzip -t / python zipfile),失败自动重试续传
set -e
BASE="$(cd "$(dirname "$0")/.." && pwd)"
RAW="$BASE/data/raw"
mkdir -p "$RAW"; cd "$RAW"

gz_ok()  { gzip -t "$1" 2>/dev/null; }
zip_ok() { python3 -c "import zipfile,sys; zipfile.ZipFile(sys.argv[1]).testzip()" "$1" >/dev/null 2>&1; }

dl() {
  url=$1; out=$2; verify=$3
  if [ -f "$out" ] && eval "$verify"; then echo "✓ $out 已存在"; return 0; fi
  for i in $(seq 1 15); do
    curl -fSL -C - --retry 2 -m 1800 -o "$out" "$url" 2>/dev/null || true
    if eval "$verify"; then echo "✓ $out ($(du -h "$out" | cut -f1))"; return 0; fi
    sleep 5
  done
  echo "✗ $out 下载失败(网络受限时可改用分段下载,见 README)"; return 1
}

GEO=https://ftp.ncbi.nlm.nih.gov/geo
dl "$GEO/series/GSE17nnn/GSE17649/matrix/GSE17649_series_matrix.txt.gz" \
   GSE17649_series_matrix.txt.gz "gz_ok GSE17649_series_matrix.txt.gz"
dl "$GEO/series/GSE111nnn/GSE111828/suppl/GSE111828_results_20180108_paracetamol_time_series.xlsx" \
   GSE111828_results.xlsx "zip_ok GSE111828_results.xlsx"
dl "$GEO/series/GSE120nnn/GSE120652/matrix/GSE120652_series_matrix.txt.gz" \
   GSE120652_series_matrix.txt.gz "gz_ok GSE120652_series_matrix.txt.gz"
dl "$GEO/platforms/GPL1nnn/GPL1261/annot/GPL1261.annot.gz" GPL1261.annot.gz "gz_ok GPL1261.annot.gz"
dl "$GEO/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz" GPL6244.annot.gz "gz_ok GPL6244.annot.gz"

# 氧化应激基因集 GO:0006979(EBI GOA 直接注释)
for sp in mouse human; do
  ucase=$(echo $sp | tr a-z A-Z)
  dl "https://ftp.ebi.ac.uk/pub/databases/GO/goa/${ucase}/goa_${sp}.gaf.gz" \
     goa_${sp}.gaf.gz "gz_ok goa_${sp}.gaf.gz"
  zcat goa_${sp}.gaf.gz 2>/dev/null | awk -F"\t" '$0 !~ /^!/ && $5=="GO:0006979" {print $3"\t"$2}' \
    | sort -u > /tmp/os_${sp}.tsv
  (echo "gene_symbol,entrez_id"; sed "s/\t/,/g" /tmp/os_${sp}.tsv) > "$BASE/data/os_genes_${sp}.csv"
  echo "✓ os_genes_${sp}.csv: $(($(wc -l < "$BASE/data/os_genes_${sp}.csv") - 1)) 个基因"
done
echo "[01] 数据下载完成"
