#!/usr/bin/env Rscript
# CD-GPT / PT 项目: Bioconductor + CRAN 依赖安装 (幂等/可续跑)
# 运行: Rscript 00_install_bioc.R  (R_LIBS_USER 已在 ~/.Renviron 指向 ~/Rlibs)
Sys.setenv(R_LIBS_USER = Sys.getenv("R_LIBS_USER", unset = paste0(Sys.getenv("HOME"), "/Rlibs")))
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
cat("User lib =", .libPaths()[1], "\n")

have <- function(p) requireNamespace(p, quietly = TRUE)
need <- function(p) !have(p)

## ---- 0. CRAN 基础依赖 ----
cran_pkgs <- c("tidyverse","data.table","RColorBrewer","corrplot","circlize","ComplexHeatmap",
               "cowplot","patchwork","viridis","caTools","survival","survminer","glmnet",
               "randomForest","e1071","caret","pROC","factoextra",
               "pheatmap","reshape2","stringr","openxlsx","writexl","devtools","remotes")
if (any(sapply(cran_pkgs, need))) {
  install.packages(cran_pkgs[sapply(cran_pkgs, need)],
                   repos = "https://cloud.r-project.org", Ncpus = 8, quiet = TRUE)
}

## ---- 1. Bioconductor 管理器 ----
if (need("BiocManager")) install.packages("BiocManager", repos = "https://cloud.r-project.org", Ncpus = 8)
suppressMessages(library(BiocManager, quietly = TRUE))
ok_bioc <- tryCatch({
  BiocManager::install(version = "3.18", ask = FALSE, update = FALSE)
  TRUE
}, error = function(e) { cat("  !! BiocManager version set warning:", conditionMessage(e), "\n"); FALSE })
if (!ok_bioc) {
  cat("  (continuing with current BiocManager repositories)\n")
}

## ---- 2. 核心生信包 (逐包 tryCatch, 失败不中断) ----
bioc_core <- c("GEOquery","limma","edgeR","DESeq2","WGCNA","ConsensusClusterPlus",
               "minfi","DMRcate","impute","sva","batcheffect",
               "org.Hs.eg.db","GSVA","GSEABase")
install_one <- function(p) {
  if (have(p)) { cat(sprintf("  [skip]   %-32s already installed\n", p)); return(TRUE) }
  cat(sprintf("  [install]%-32s ...\n", p))
  ok <- tryCatch({ BiocManager::install(p, ask = FALSE, update = FALSE, Ncpus = 8); TRUE },
                 error = function(e){ cat("    !! error:", conditionMessage(e), "\n"); FALSE })
  if (ok && have(p)) { cat(sprintf("  [done]   %-32s OK\n", p)); TRUE }
  else { cat(sprintf("  [FAIL]   %-32s NOT installed\n", p)); FALSE }
}
for (p in bioc_core) install_one(p)

## ---- 3. EPIC 芯片注释/平台包 (大, 单独处理) ----
anno_pkgs <- c("IlluminaHumanMethylationEPICmanifest",
               "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
               "IlluminaHumanMethylation450kmanifest",
               "IlluminaHumanMethylation450kanno.ilmn12.hg19",
               "TxDb.Hsapiens.UCSC.hg19.knownGene")
for (p in anno_pkgs) install_one(p)

## ---- 4. 汇总 ----
core_needed <- c("minfi","ConsensusClusterPlus",
                 "org.Hs.eg.db","GEOquery","limma","WGCNA","DMRcate","impute",
                 "IlluminaHumanMethylationEPICmanifest",
                 "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                 "IlluminaHumanMethylation450kmanifest",
                 "IlluminaHumanMethylation450kanno.ilmn12.hg19")
cat("\n===== INSTALL SUMMARY =====\n")
cat("Total packages in user lib:", length(list.files(.libPaths()[1])), "\n")
missing0 <- core_needed[sapply(core_needed, need)]
if (length(missing0) == 0) {
  cat("ALL CORE PACKAGES PRESENT\n")
} else {
  cat("STILL MISSING:", paste(missing0, collapse = ", "), "\n")
}

## ---- 5. 缺失包重试 (最多 3 轮) ----
attempt <- 0
while (attempt < 3) {
  missing_now <- core_needed[sapply(core_needed, need)]
  if (length(missing_now) == 0) break
  attempt <- attempt + 1
  cat(sprintf("\n--- retry round %d for: %s ---\n", attempt, paste(missing_now, collapse = ", ")))
  for (p in missing_now) install_one(p)
}
missing_final <- core_needed[sapply(core_needed, need)]
if (length(missing_final) == 0) {
  cat("\nALL CORE PACKAGES PRESENT - writing install.done\n")
  writeLines(c(format(Sys.time()), "all core packages present"),
             file.path(.libPaths()[1], "..", "install.done"))
} else {
  cat("\nWARNING - still missing after retries:", paste(missing_final, collapse = ", "), "\n")
  cat("Please re-run this script.\n")
}
