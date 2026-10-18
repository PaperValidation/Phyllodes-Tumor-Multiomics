#!/usr/bin/env Rscript
# =====================================================================
# 00_download_geoTumor.R
# 乳腺叶状肿瘤 (Breast Phyllodes Tumor) 项目 — GEO 数据获取与样本表重建
#
# 用法 (在 WSL 上):
#   Rscript 00_download_geoTumor.R "/mnt/e/Medical Data/Dermatology/Phyllodes tumor"
# 或不带参数, 读取环境变量 DATASET_PATH; 再不行用默认路径。
#
# 依赖 R 包: GEOquery, Biobase, data.table
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install(c("GEOquery","Biobase"))
#   install.packages("data.table")
#
# 说明:
#   - getGSE 默认只下载 family SOFT + 系列矩阵(表达式/beta 值), 不会自动拉 IDAT,
#     对 96 例 EPIC 的 GSE179458 而言体量与速度都可控。
#   - 若后续需要 IDAT 级 minfi 处理, 再单独用 GEOquery::getGEOSuppFiles() 拉补充文件。
# =====================================================================

args <- commandArgs(trailingOnly = TRUE)
data_root <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  Sys.getenv("DATASET_PATH", "/mnt/e/Medical Data/Dermatology/Phyllodes tumor")
}

## 四个核心 GEO 系列 (均开放)
gse_ids <- c("GSE179458", "GSE231574", "GSE235267", "GSE78071")

## 各系列在下游流程中的角色
gse_role <- c(
  GSE179458 = "methylation_core",    # 96 例 EPIC 甲基化 + CNV (benign/borderline/malignant)
  GSE231574 = "methylation_second",  # 33 例 EPIC 甲基化
  GSE235267 = "lncRNA_mRNA",         # 16 例 lncRNA + mRNA 芯片
  GSE78071  = "expression_FA_PT"     # 48 例 FA/PT 表达芯片
)

cat(sprintf("[INFO] data_root = %s\n", data_root))
if (!dir.exists(data_root)) {
  dir.create(data_root, recursive = TRUE)
  cat(sprintf("[INFO] created data_root\n"))
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(data.table)
})

download_one <- function(gse, root) {
  dest <- file.path(root, "raw", gse)
  if (!dir.exists(dest)) dir.create(dest, recursive = TRUE)
  cat(sprintf("\n[DOWNLOAD] %s (%s) -> %s\n",
              gse, gse_role[[gse]], dest))

  gse_obj <- getGEO(gse, destdir = dest, AnnotGPL = FALSE)
  if (is.list(gse_obj) && !inherits(gse_obj, "GSE")) gse_obj <- gse_obj[[1]]

  ## 保存完整 GSE 对象
  saveRDS(gse_obj, file.path(dest, sprintf("%s.rds", gse)))

  ## 提取并保存表型 (pData) — 下游据此重建 benign/borderline/malignant/FA 标签
  pd <- as.data.frame(pData(phenoData(gse_obj)))
  fwrite(pd, file.path(dest, sprintf("%s_pdata.csv", gse)))
  cat(sprintf("[INFO] pData 列名: %s\n",
              paste(colnames(pd), collapse = ", ")))

  ## 保存各平台表达/beta 矩阵
  elems <- names(gse_obj@gset)
  cat(sprintf("[INFO] %s 元素: %s\n", gse, paste(elems, collapse = ", ")))
  for (nm in elems) {
    m <- tryCatch(exprs(gse_obj@gset[[nm]]), error = function(e) NULL)
    if (!is.null(m) && is.matrix(m)) {
      write.csv(m, file.path(dest, sprintf("%s_%s_expr.csv", gse, nm)))
      cat(sprintf("[INFO] 已存矩阵 %s (%d features x %d samples)\n",
                  nm, nrow(m), ncol(m)))
    }
  }
  invisible(gse)
}

for (g in gse_ids) download_one(g, data_root)

cat("\n[DONE] 四个系列已下载至", file.path(data_root, "raw"), "\n")
cat("[NEXT] 检查各 *_pdata.csv, 按 characteristics_ch1 重建样本层级标签,\n")
cat("       再进入 01_methylation_pipeline.R (GSE179458+GSE231574 共识聚类 + DMP/DMR)。\n")
