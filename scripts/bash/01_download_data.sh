#!/usr/bin/env bash
# 断点续跑式 GEO 下载脚本 (checkpoint-based)
# 用法: bash 01_download_data.sh
# 每个文件下载完成写 .done 标记; 重跑自动跳过已完成文件, 未完成的用 curl -C - 续传
set +e
# 数据根: 从 .env 解析 (DATASET_PATH 为唯一真相源)
ENV_FILE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/.env"
_DP="$(grep -E '^DATASET_PATH=' "$ENV_FILE" 2>/dev/null | head -1 | cut -d= -f2- | tr -d "\"'")"
BASE="${_DP//\\//}"
case "$BASE" in
  [A-Za-z]:/*) BASE="/mnt/$(printf '%s' "${BASE:0:1}" | tr 'A-Z' 'a-z')/${BASE:2}" ;;
esac
RAW="$BASE/raw"
PROC="$BASE/data/processed"
LOG="$PROC/download_data.log"
mkdir -p "$RAW" "$PROC"
echo "===== download start $(date) =====" >> "$LOG"

download_one () {
  # $1 = URL, $2 = 目标路径
  local url="$1" dest="$2"
  local mark="${dest}.done"
  if [ -f "$mark" ]; then
    echo "  [skip, done] $(basename "$dest")" | tee -a "$LOG"
    return 0
  fi
  echo ">> $(date +%H:%M:%S) downloading: $(basename "$dest")" | tee -a "$LOG"
  # 续传最多重试8次; 完成后写 .done
  for i in $(seq 1 8); do
    curl -sL --retry 3 -C - -o "$dest" "$url" >> "$LOG" 2>&1
    if [ -s "$dest" ]; then
      # 简单完整性校验: 对压缩包尝试列出内容
      case "$dest" in
        *.tar) tar -tf "$dest" >/dev/null 2>&1 && { touch "$mark"; echo "   [ok] $(basename "$dest") ($(du -h "$dest"|cut -f1))" | tee -a "$LOG"; return 0; } ;;
        *.gz)  gzip -t "$dest" >/dev/null 2>&1 && { touch "$mark"; echo "   [ok] $(basename "$dest") ($(du -h "$dest"|cut -f1))" | tee -a "$LOG"; return 0; } ;;
        *)     touch "$mark"; echo "   [ok] $(basename "$dest") ($(du -h "$dest"|cut -f1))" | tee -a "$LOG"; return 0 ;;
      esac
    fi
    echo "   retry $i ..." | tee -a "$LOG"
    sleep 3
  done
  echo "   [FAILED] $(basename "$dest")" | tee -a "$LOG"
  return 1
}

# ---------- GSE231574 (33例 EPIC 已处理 beta, 先校验跳过) ----------
D="$RAW/GSE231574"; mkdir -p "$D"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231574/suppl/GSE231574_beta_detp.csv.gz" "$D/GSE231574_beta_detp.csv.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231574/suppl/GSE231574_samplesPhyllodes.csv.gz" "$D/GSE231574_samplesPhyllodes.csv.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231574/suppl/GSE231574_raw_detp.csv.gz" "$D/GSE231574_raw_detp.csv.gz"

# ---------- GSE235267 (16例 lncRNA+mRNA 芯片, 小优先) ----------
D="$RAW/GSE235267"; mkdir -p "$D"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235267/matrix/GSE235267_series_matrix.txt.gz" "$D/GSE235267_series_matrix.txt.gz"

# ---------- GSE78071 (48例 FA/PT 表达芯片, 小优先) ----------
D="$RAW/GSE78071"; mkdir -p "$D"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78071/matrix/GSE78071_series_matrix.txt.gz" "$D/GSE78071_series_matrix.txt.gz"

# ---------- GSE179458 (96例 EPIC IDAT 原始, 1.57GB, 最后续传) ----------
D="$RAW/GSE179458"; mkdir -p "$D"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179458/suppl/GSE179458_RAW.tar" "$D/GSE179458_RAW.tar"
# 03 解析表型真值(!Sample_title) 需要 series matrix
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179458/matrix/GSE179458_series_matrix.txt.gz" "$D/GSE179458_series_matrix.txt.gz"

echo "===== download end $(date) =====" >> "$LOG"
echo "DOWNLOAD_ALL_DONE" >> "$LOG"
