#!/bin/bash
# 00_install_deps.sh — 依赖自检与补装
# 校验 rna 环境;org.Mm/org.Hs.eg.db 为 Bioconductor 大数据包,conda 通道
# 仅含 stub,故在此用分段下载源码安装(单连接下载在网络受限环境必失败)
set -e
BASE="$(cd "$(dirname "$0")/.." && pwd)"
source ~/miniforge3/bin/activate rna 2>/dev/null || {
  echo "未找到 rna 环境,正在创建..."
  mamba env create -f "$BASE/env/rna.yml"
  source ~/miniforge3/bin/activate rna
}
LIB=$(Rscript -e "cat(.libPaths()[1])")

# 分段并行下载:精确到每段字节数的校验(≥ 会导致错位,必须 ==)
seg_dl() {
  url=$1; out=$2; n=${3:-12}
  cl=$(curl -sIL -m 30 "$url" | grep -i "^content-length" | tail -1 | tr -dc "0-9")
  [ -z "$cl" ] && return 1
  chunk=$(( (cl + n - 1) / n ))
  for ((s=0; s<n; s++)); do
    ( st=$(( s * chunk )); en=$(( st + chunk - 1 )); [ $en -ge $cl ] && en=$(( cl - 1 ))
      want=$(( en - st + 1 ))
      for t in $(seq 1 25); do
        have=0; [ -f "$out.p$s" ] && have=$(stat -c %s "$out.p$s")
        [ "$have" -eq "$want" ] && break
        [ "$have" -gt "$want" ] && { rm -f "$out.p$s"; have=0; }
        curl -fsSL -m 600 -r $(( st + have ))-${en} -o - "$url" >> "$out.p$s" 2>/dev/null
        sleep 2
      done ) &
  done
  wait
  total=0
  for ((s=0; s<n; s++)); do [ -f "$out.p$s" ] && total=$(( total + $(stat -c %s "$out.p$s") )); done
  [ "$total" -ne "$cl" ] && return 1
  for ((s=0; s<n; s++)); do cat "$out.p$s"; done > "$out"
  rm -f "$out".p*
  return 0
}

install_db() {
  Rscript -e "cat(requireNamespace('$1', quietly=TRUE))" | grep -q TRUE && { echo "✓ $1 已安装"; return 0; }
  echo "安装 $1(分段下载 ~90MB,约 3-10 分钟)..."
  for u in "https://mghp.osn.xsede.org/bir190004-bucket01/archive.bioconductor.org/packages/3.22/data/annotation/src/contrib/${1}_${2}.tar.gz" \
           "https://bioconductor.org/packages/3.22/data/annotation/src/contrib/${1}_${2}.tar.gz"; do
    if seg_dl "$u" "/tmp/${1}.tar.gz" 16; then
      if R CMD INSTALL -l "$LIB" "/tmp/${1}.tar.gz" >/dev/null 2>&1; then
        echo "✓ $1 安装成功"; rm -f "/tmp/${1}.tar.gz"; return 0
      fi
    fi
  done
  echo "✗ $1 安装失败"; return 1
}

install_db org.Mm.eg.db 3.22.0
install_db org.Hs.eg.db 3.22.0
echo "[00] 依赖就绪"
