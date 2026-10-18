#!/usr/bin/env Rscript
# ============================================================
# 02_meth_gse231574.R — 甲基化主队列分析 (GSE231574, 33例 EPIC)
# 输入: raw/GSE231574/GSE231574_beta_detp.csv.gz (宽表: beta_<SID>, p_det_p_det_<SID>)
#       raw/GSE231574/GSE231574_samplesPhyllodes.csv.gz (barcode, Sample ID, SGH Score, New Classification grade)
# 输出: data/processed/gse231574_*.rds + results/tables + results/figures
# 流程: 加载β->表型(Sample ID对齐)->探针过滤->共识聚类->DMP(limma)->DMR(DMRcate)->GO
# ============================================================
suppressMessages({
  library(minfi); library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(ConsensusClusterPlus); library(limma); library(DMRcate)
  library(org.Hs.eg.db); library(AnnotationDbi); library(ggplot2)
})
# 项目路径唯一真相源 (对齐 scripts/_env.py)
.cmdArgs <- commandArgs(trailingOnly = FALSE)
.fileArg <- sub("^--file=", "", .cmdArgs[grepl("^--file=", .cmdArgs)])
.sdir <- if (length(.fileArg)) normalizePath(dirname(.fileArg[1]), winslash="/") else normalizePath(getwd(), winslash="/")
source(file.path(.sdir, "_env.R"))
source(file.path(.sdir, "enrich_helpers.R"))

DS  <- DATASET_PATH
RAW <- file.path(DS, "raw/GSE231574")
PRC <- file.path(DS, "data/processed"); if (!dir.exists(PRC)) dir.create(PRC, recursive = TRUE)
RES <- file.path(DS, "results"); if (!dir.exists(RES)) dir.create(RES, recursive = TRUE)
TAB <- file.path(RES, "tables"); FIG <- file.path(RES, "figures")
dir.create(TAB, recursive = TRUE); dir.create(FIG, recursive = TRUE)

cat(">> [1] 加载 β 宽表 (beta_<SID> / p_det_p_det_<SID> 交错)\n")
mat <- read.csv(gzfile(file.path(RAW, "GSE231574_beta_detp.csv.gz")),
                row.names = 1, check.names = FALSE)
beta_cols <- grep("^beta_", colnames(mat), value = TRUE)
detp_cols <- grep("^p_det_p_det_", colnames(mat), value = TRUE)
beta <- as.matrix(mat[, beta_cols, drop = FALSE])
detp <- as.matrix(mat[, detp_cols, drop = FALSE])
sid_beta <- sub("^beta_", "", beta_cols)          # "4018"
sid_detp <- sub("^p_det_p_det_", "", detp_cols)   # "4018"
stopifnot(all.equal(sid_beta, sid_detp))          # beta 与 detp 严格成对
colnames(beta) <- sid_beta; colnames(detp) <- sid_detp
cat(sprintf("   β: %d probes x %d samples\n", nrow(beta), ncol(beta)))

cat(">> [2] 表型整理 (按 Sample ID 对齐)\n")
phen <- read.csv(gzfile(file.path(RAW, "GSE231574_samplesPhyllodes.csv.gz")),
                 stringsAsFactors = FALSE, check.names = FALSE)
phen$Sample.ID <- as.character(phen$`Sample ID`)
common <- intersect(colnames(beta), phen$Sample.ID)
cat(sprintf("   匹配到 %d/%d 样本 (phen 共 %d 行)\n", length(common), ncol(beta), nrow(phen)))
beta <- beta[, common, drop = FALSE]; detp <- detp[, common, drop = FALSE]
phen <- phen[match(common, phen$Sample.ID), ]
phen$group_raw <- phen$`New Classification grade`
phen$group3 <- ifelse(phen$group_raw=="Fibroadenoma","FA",
              ifelse(phen$group_raw=="Metaplastic","Metaplastic", as.character(phen$group_raw)))
sev_map <- c("Benign"=1, "Borderline"=2, "Malignant"=3)
phen$severity <- sev_map[phen$group3]
phen$is_PT <- phen$group_raw %in% c("Benign","Borderline","Malignant")
saveRDS(phen, file.path(PRC, "gse231574_pheno.rds"))
cat("   分组计数:\n"); print(table(phen$group3))

cat(">> [3] 探针过滤\n")
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep <- intersect(rownames(beta), rownames(anno))
keep <- keep[anno[keep, "chr"] %in% paste0("chr", 1:22)]
for (sc in c("SNP_10","Probe_rs")) if (sc %in% colnames(anno)) keep <- keep[is.na(anno[keep, sc])]
xr_url <- "https://www.txwikinger.me.uk/crosshyb/CrossHyb_EPIC.zip"
xr_file <- file.path(PRC, "CrossHyb_EPIC.txt")
tryCatch({
  if (!file.exists(xr_file)) download.file(xr_url, xr_file, mode="wb", quiet=TRUE)
  xr <- read.table(xr_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  xr_probes <- unique(c(as.character(xr$TargetID), as.character(xr$CrossHybID)))
  keep <- setdiff(keep, xr_probes); cat(sprintf("   去除交叉反应探针, 余 %d\n", length(keep)))
}, error = function(e) cat("   (交叉反应列表下载失败, 跳过该步)\n"))
keep <- keep[rowSums(detp[keep, , drop=FALSE] > 0.01) == 0]
keep <- keep[apply(beta[keep, , drop=FALSE], 1, function(r) all(r>=0 & r<=1, na.rm=TRUE))]
beta_f <- beta[keep, ]; beta_f[is.na(beta_f)] <- 0.5
anno_f <- anno[keep, ]
cat(sprintf("   最终保留探针: %d (起始 %d)\n", nrow(beta_f), nrow(beta)))
saveRDS(beta_f, file.path(PRC, "gse231574_beta_filtered.rds"))
saveRDS(anno_f, file.path(PRC, "gse231574_anno_filtered.rds"))

cat(">> [4] M值转换 (聚类/DMP 用)\n")
Mval <- log2(beta_f/(1-beta_f+1e-6))
saveRDS(Mval, file.path(PRC, "gse231574_Mval_filtered.rds"))

cat(">> [5] 共识聚类 (PT 进展谱系内, K=2 高风险/低风险二分)\n")
set.seed(20240728)
pt_idx <- which(phen$is_PT)
M_pt <- Mval[, pt_idx, drop=FALSE]
# 取 PT 样本间变异最大的前 5000 探针聚类 (50万探针聚类无意义且爆内存)
mv <- apply(M_pt, 1, mad, na.rm=TRUE); topf <- head(order(-mv), 5000)
X <- t(M_pt[topf, , drop=FALSE])            # n_PT 样本 × 5000 探针
d <- dist(X, method="euclidean")
hc <- hclust(d, method="complete")
k <- 2
clust <- cutree(hc, k=k)
phen$cluster_pt <- NA; phen$cluster_pt[pt_idx] <- clust
saveRDS(phen, file.path(PRC, "gse231574_pheno.rds"))
cat(sprintf("   K=%d 聚类分布: %s (基于 %d 高变异探针, hclust/euclidean)\n", k, paste(table(clust), collapse="/"), length(topf)))
# 树状图 + 热图 (稿件图)
png(file.path(FIG, "gse231574_cluster_heatmap.png"), width=1100, height=900, res=150)
heatmap(t(X), scale="row", col=colorRampPalette(c("blue","white","red"))(50),
        labRow=FALSE, labCol=phen$group3[pt_idx],
        main="GSE231574 PT subtypes (top5000 MAD probes)")
dev.off()

cat(">> [6] DMP 分析 (limma: 进展梯度 + 组间两两, 仅 PT)\n")
grp <- factor(phen$group3[pt_idx], levels=c("Benign","Borderline","Malignant"))
design <- model.matrix(~0+grp); colnames(design) <- levels(grp)
fit <- lmFit(M_pt, design)
contrasts <- makeContrasts(
  Mal_vs_Ben = Malignant-Benign,
  Bor_vs_Ben = Borderline-Benign,
  Mal_vs_Bor = Malignant-Borderline,
  Trend      = (Borderline-Benign)+(Malignant-Borderline),
  levels=design)
fit2 <- contrasts.fit(fit, contrasts); fit2 <- eBayes(fit2)
for (cn in colnames(contrasts)) {
  tt <- topTable(fit2, coef=cn, number=Inf, p.value=1)
  tt$probe <- rownames(tt)
  tt$gene <- anno_f[tt$probe, "UCSC_RefGene_Name"]
  tt$region <- anno_f[tt$probe, "Regulatory_Feature_Group"]
  write.csv(tt, file.path(TAB, sprintf("gse231574_DMP_%s.csv", cn)), row.names=FALSE)
  cat(sprintf("   DMP %-12s: %d 显著 (FDR<0.05)\n", cn, sum(tt$adj.P.Val<0.05)))
}

cat(">> [7] DMR 分析 (DMRcate, 恶性进展趋势 severity)\n")
myAnnotation <- anno_f[, c("chr","pos","strand")]; colnames(myAnnotation) <- c("chr","pos","strand")
rownames(myAnnotation) <- rownames(anno_f)
sev_pt <- phen$severity[pt_idx]
design.dmr <- cbind("(Intercept)"=1, severity=sev_pt)
myBetas <- Mval[, pt_idx, drop=FALSE]   # cpg.annotate what="M" 需 M 值 (非 β)
dmr <- cpg.annotate("array", myBetas, what="M", arraytype="EPIC",
                    analysis.type="differential", design=design.dmr, coef=2, anno=myAnnotation)
nsig <- if ("fdr" %in% colnames(dmr)) sum(dmr$fdr < 0.05, na.rm=TRUE) else 0
cat(sprintf("   cpg.annotate 显著CpG (FDR<0.05): %d\n", nsig))
dmr.dt <- data.frame()
if (nsig > 0) {
  dmr.res <- dmrcate(dmr, lambda=1000, C=2)
  dmr.dt <- data.frame(extractRanges(dmr.res, genome="hg19"))
  dmr.dt <- dmr.dt[dmr.dt$n.sig>0, ]
  cat(sprintf("   DMR (FDR): %d 个区域 (含显著CpG)\n", nrow(dmr.dt)))
} else {
  # FDR 无显著: 用未校正 p<0.05 提取候选 DMR (探索性, 供稿件讨论)
  dmr.res <- tryCatch(dmrcate(dmr, lambda=1000, C=2, pcutoff=0.001), error=function(e) NULL)
  if (!is.null(dmr.res)) {
    dmr.dt <- data.frame(extractRanges(dmr.res, genome="hg19"))
    dmr.dt <- dmr.dt[dmr.dt$n.sig>0, ]
    if (nrow(dmr.dt) > 0) dmr.dt$threshold <- "raw_p<0.001_exploratory"
    cat(sprintf("   DMR (raw p<0.001, 探索性): %d 个区域\n", nrow(dmr.dt)))
  } else {
    cat("   (FDR 与 uncorrected 均无显著 DMR, 记录空表)\n")
  }
}
write.csv(dmr.dt, file.path(TAB, "gse231574_DMR_Mal_vs_Ben.csv"), row.names=FALSE)

cat(">> [8] 功能富集 (Malignant vs Benign DMP)\n")
dmp_mb <- read.csv(file.path(TAB, "gse231574_DMP_Mal_vs_Ben.csv"))
sig <- dmp_mb[dmp_mb$adj.P.Val<0.05, ]
genes <- unique(unlist(strsplit(as.character(sig$gene), ";")))
genes <- genes[genes!="" & !is.na(genes)]
if (length(genes) > 0) {
  eg <- AnnotationDbi::select(org.Hs.eg.db, keys=genes, columns="ENTREZID", keytype="SYMBOL")
  eg <- eg[!is.na(eg$ENTREZID), ]
  go <- enrichGO_manual(eg$ENTREZID, ont="ALL",
                        pAdjustMethod="BH", pvalueCutoff=0.05)
  write.csv(as.data.frame(go), file.path(TAB, "gse231574_DMP_Mal_vs_Ben_GO.csv"), row.names=FALSE)
  cat(sprintf("   GO 富集条目: %d (输入基因 %d)\n", ifelse(is.null(go),0,nrow(go)), length(genes)))
} else cat("   (无显著 DMP 基因, 跳过富集)\n")

cat(">> DONE: 02_meth_gse231574.R\n")
