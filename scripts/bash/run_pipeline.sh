#!/usr/bin/env bash
# ============================================================
# PT 项目全流程编排 (run_pipeline.sh)
# 作用: 确保 Bioconductor 包就绪 -> 顺序跑 02-06 (断点续跑)
# 运行: wsl -e bash /mnt/c/.../scripts/bash/run_pipeline.sh
# 日志: data/processed/pipeline_run.log (Rscript 输出) + harness stdout (进度)
# ============================================================
set -u
# 项目根: 由脚本位置推导 (不再硬编码 C: 路径)
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# 数据根: 从 .env 解析 (DATASET_PATH 为唯一真相源)
ENV_FILE="$BASE/.env"
_DP="$(grep -E '^DATASET_PATH=' "$ENV_FILE" 2>/dev/null | head -1 | cut -d= -f2- | tr -d "\"'")"
DS="${_DP//\\//}"
case "$DS" in
  [A-Za-z]:/*) DS="/mnt/$(printf '%s' "${DS:0:1}" | tr 'A-Z' 'a-z')/${DS:2}" ;;
esac
RAW="$DS/raw"; PRC="$DS/data/processed"; RES="$DS/results"
SCR="$BASE/scripts/R"; TAB="$RES/tables"; FIG="$RES/figures"
LOG="$PRC/pipeline_run.log"
export R_LIBS_USER="/home/tiduslrg/Rlibs"
mkdir -p "$TAB" "$FIG"
ts(){ date '+%Y-%m-%d %H:%M:%S'; }
echo "[$(ts)] ===== PT pipeline orchestrator START ====="
: > "$LOG"   # 每次启动清空日志, 便于监控本论运行

# ---------- 0. 确保所需 Bioconductor 包真正就绪 ----------
REQ=(minfi ConsensusClusterPlus org.Hs.eg.db \
      GEOquery limma WGCNA DMRcate impute IlluminaHumanMethylationEPICanno.ilm10b4.hg19 \
      IlluminaHumanMethylation450kmanifest IlluminaHumanMethylation450kanno.ilmn12.hg19)
MAXWAIT=21600   # 6h 上限
WAITED=0
pkg_ok(){
  for p in "${REQ[@]}"; do
    if ! Rscript -e "cat(any(dir.exists(file.path(.libPaths(), '$p'))))" 2>/dev/null | grep -q TRUE; then
      echo "$p"
      return 1
    fi
  done
  return 0
}
install_done_marker="/home/tiduslrg/install.done"
while true; do
  miss_list="$(pkg_ok)" || true
  if [ -z "$miss_list" ]; then
    echo "[$(ts)] all required packages present"; break
  fi
  echo "[$(ts)] missing packages: $miss_list"
  if [ "$WAITED" -ge "$MAXWAIT" ]; then
    echo "[$(ts)] TIMEOUT waiting for packages"; exit 2
  fi
  # 若安装进程在跑, 等它; 否则主动重跑安装
  if pgrep -f "00_install_bioc.R" >/dev/null 2>&1; then
    echo "[$(ts)] install process running, wait 60s"; sleep 60; WAITED=$((WAITED+60)); continue
  fi
  # 无安装进程在跑: 删除可能存在的失效 install.done, 启动安装(foreground, 由本脚本托管)
  if [ -f "$install_done_marker" ]; then
    echo "[$(ts)] install.done stale (pkg missing) -> remove"; rm -f "$install_done_marker"
  fi
  echo "[$(ts)] launching install (foreground)..."
  Rscript "$SCR/00_install_bioc.R" >> "$PRC/install_bioc_run.log" 2>&1
  echo "[$(ts)] install run finished, re-checking packages"; sleep 5
done

# ---------- 1. 顺序跑分析 (断点续跑) ----------
run_step(){
  local id="$1" script="$2"
  local donef="$PRC/${id}.done"
  if [ -f "$donef" ]; then echo "[$(ts)] $id: already done (skip)"; return 0; fi
  echo "[$(ts)] >>> $id: Rscript $script"
  if Rscript "$SCR/$script" >> "$LOG" 2>&1; then
    touch "$donef"; echo "[$(ts)] <<< $id: OK"; return 0
  else
    echo "[$(ts)] !!! $id: FAILED (see $LOG)"; return 1
  fi
}

run_step 02_meth_gse231574     02_meth_gse231574.R      || { echo "[$(ts)] ABORT after 02"; exit 1; }
run_step 03_meth_gse179458     03_meth_gse179458.R      || { echo "[$(ts)] ABORT after 03"; exit 1; }
run_step 04_lncrna_gse235267   04_lncrna_gse235267.R    || { echo "[$(ts)] ABORT after 04"; exit 1; }
run_step 05_expr_gse78071      05_expr_gse78071.R       || { echo "[$(ts)] ABORT after 05"; exit 1; }
run_step 06_drug_repositioning 06_drug_repositioning.R  || { echo "[$(ts)] ABORT after 06"; exit 1; }

echo "[$(ts)] ===== PT pipeline COMPLETE ====="
