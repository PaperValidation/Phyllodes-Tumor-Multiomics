#!/usr/bin/env Rscript
# ============================================================
# 03_meth_gse179458.R — 甲基化验证队列 (GSE179458, 96例 EPIC IDAT)
# 输入: raw/GSE179458/*.idat.gz (已解tar) + series_matrix.txt.gz (表型真值)
# 流程: 解压IDAT -> minfi读 -> funnorm -> β/M -> 表型解析 -> 探针过滤
#       -> 共识聚类 -> DMP (恶性vs良性PT) -> 与主队列交叉验证
# ============================================================
suppressMessages({
  library(minfi); library(IlluminaHumanMethylation450kmanifest)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(ConsensusClusterPlus); library(limma); library(DMRcate)
})

# 项目路径唯一真相源 (对齐 scripts/_env.py)
.cmdArgs <- commandArgs(trailingOnly = FALSE)
.fileArg <- sub("^--file=", "", .cmdArgs[grepl("^--file=", .cmdArgs)])
.sdir <- if (length(.fileArg)) normalizePath(dirname(.fileArg[1]), winslash="/") else normalizePath(getwd(), winslash="/")
source(file.path(.sdir, "_env.R"))

DS  <- DATASET_PATH
RAW <- file.path(DS, "raw/GSE179458")
PRC <- file.path(DS, "data/processed"); dir.create(PRC, recursive=TRUE)
RES <- file.path(DS, "results"); TAB <- file.path(RES, "tables")
dir.create(TAB, recursive=TRUE)

cat(">> [1] 解压 RAW.tar -> IDAT.gz -> IDAT\n")
tarf <- file.path(RAW, "GSE179458_RAW.tar")
if (file.exists(tarf)) {
  cat("   解压 GSE179458_RAW.tar ...\n")
  untar(tarf, exdir = RAW)   # 内置 untar 自动处理含空格路径
}
gz <- list.files(RAW, pattern="\\.idat\\.gz$", full.names=TRUE, recursive=TRUE)
if (length(gz) > 0) {
  cat(sprintf("   解压 %d 个 idat.gz ...\n", length(gz)))
  invisible(sapply(gz, function(f) system(paste0("gunzip -f '", f, "'"), ignore.stderr=TRUE)))
}

cat(">> [2] minfi 读取并归一化\n")
idat_dir <- RAW
rg <- read.metharray.exp(base = idat_dir, recursive = FALSE,
                         verbose = FALSE, force = TRUE)   # 450k(GPL21145); force 兼容个别异常idat
cat(sprintf("   RGChannelSet: %d probes x %d samples\n", nrow(rg), ncol(rg)))
qnt <- preprocessQuantile(rg, fixOutliers=TRUE, verbose=FALSE)
beta <- getBeta(qnt); Mval <- getM(qnt)
# minfi 列名已是完整样本ID (如 203293640041_R07C01), 不做截断

cat(">> [3] 从 series_matrix 解析表型真值\n")
mtxt <- gzfile(file.path(RAW, "GSE179458_series_matrix.txt.gz"))
lines <- readLines(mtxt, warn=FALSE)
# !Sample_geo_accession 行给出 GSM 顺序
acc_line <- lines[grep("^!Sample_geo_accession", lines)[1]]
ttl_line <- lines[grep("^!Sample_title", lines)[1]]
accs <- gsub('"', '', unlist(strsplit(acc_line, "\t"))[-1])
ttls <- gsub('"', '', unlist(strsplit(ttl_line, "\t"))[-1])
names(ttls) <- accs
# minfi 列名形如 203293640041_R07C01 (Sample_title 前缀)
sample_labels <- sapply(colnames(beta), function(b) {
  bb <- sub("^GSM[0-9]+_", "", b)   # idat 文件名带 GSM 前缀, 标题里无 -> 匹配前剥掉
  m <- ttls[grep(bb, ttls, fixed=TRUE)]; ifelse(length(m), m, NA)
})
phen <- data.frame(barcode = colnames(beta), raw_label = sample_labels,
                   stringsAsFactors=FALSE)
# 分类: benign/borderline/malignant PT, PT NOS, FA(BR_FAD), BR_CA, DLBCL, 其他
phen$class <- NA
phen$class[grepl("benign PT|PHYT_BEN|benign or borderline", phen$raw_label, ignore.case=TRUE)] <- "Benign"
phen$class[grepl("borderline PT|PHYT_BOR", phen$raw_label, ignore.case=TRUE)] <- "Borderline"
phen$class[grepl("malignant PT|PHYT_MAL", phen$raw_label, ignore.case=TRUE)] <- "Malignant"
phen$class[grepl("PT NOS|PHYT_NOS", phen$raw_label, ignore.case=TRUE)] <- "PT_NOS"
phen$class[grepl("BR_FAD", phen$raw_label)] <- "Fibroadenoma"
phen$class[grepl("BR_CA", phen$raw_label)] <- "BreastCarcinoma"
phen$class[grepl("DLBCL", phen$raw_label)] <- "DLBCL"
phen$class[grepl("DNADEG", phen$raw_label)] <- "DNAdeg"
phen$class[grepl("\\bAS\\b", phen$raw_label)] <- "Adenosis"
phen$class[grepl("^-$", phen$raw_label)] <- "Unknown"
phen$is_PT <- phen$class %in% c("Benign","Borderline","Malignant")
saveRDS(phen, file.path(PRC, "gse179458_pheno.rds"))
cat("   分类计数:\n"); print(table(phen$class, useNA="ifany"))

cat(">> [4] 探针过滤 (同 02 口径)\n")
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep <- intersect(rownames(beta), rownames(anno))
keep <- keep[anno[keep,"chr"] %in% paste0("chr",1:22)]
for (sc in c("SNP_10","Probe_rs")) if (sc %in% colnames(anno)) keep <- keep[is.na(anno[keep,sc])]
xr_file <- file.path(PRC, "CrossHyb_EPIC.txt")
if (file.exists(xr_file)) {
  xr <- read.table(xr_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  keep <- setdiff(keep, unique(c(as.character(xr$TargetID), as.character(xr$CrossHybID))))
}
keep <- keep[apply(beta[keep,,drop=FALSE],1,function(r) all(r>=0&r<=1,na.rm=TRUE))]
beta_f <- beta[keep,]; beta_f[is.na(beta_f)] <- 0.5
Mval_f <- log2(beta_f/(1-beta_f+1e-6))
saveRDS(beta_f, file.path(PRC, "gse179458_beta_filtered.rds"))
saveRDS(Mval_f, file.path(PRC, "gse179458_Mval_filtered.rds"))
cat(sprintf("   保留探针 %d\n", nrow(beta_f)))

cat(">> [5] 共识聚类 (PT 进展谱系内)\n")
set.seed(20240728)
pt_idx <- which(phen$is_PT)
M_pt <- Mval_f[, pt_idx, drop=FALSE]
# 取 PT 样本间变异最大的前 5000 探针聚类 (50万探针聚类无意义且爆内存)
mv <- apply(M_pt, 1, mad, na.rm=TRUE); topf <- head(order(-mv), 5000)
X <- t(M_pt[topf, , drop=FALSE])            # n_PT 样本 × 5000 探针
d <- dist(X, method="euclidean")
hc <- hclust(d, method="complete")
k <- 2
clust <- cutree(hc, k=k)
phen$cluster_pt <- NA; phen$cluster_pt[pt_idx] <- clust
saveRDS(phen, file.path(PRC, "gse179458_pheno.rds"))
cat(sprintf("   K=%d 聚类分布: %s (基于 %d 高变异探针, hclust/euclidean)\n", k, paste(table(clust), collapse="/"), length(topf)))
FIG <- file.path(RES, "figures"); dir.create(FIG, recursive=TRUE)
png(file.path(FIG, "gse179458_cluster_heatmap.png"), width=1100, height=900, res=150)
heatmap(t(X), scale="row", col=colorRampPalette(c("blue","white","red"))(50),
        labRow=FALSE, labCol=phen$class[pt_idx],
        main="GSE179458 PT subtypes (top5000 MAD probes)")
dev.off()

cat(">> [6] DMP (Malignant vs Benign PT)\n")
design <- model.matrix(~0+factor(phen$class[pt_idx]))
colnames(design) <- make.names(unique(phen$class[pt_idx]))
fit <- lmFit(M_pt, design)
cm <- makeContrasts(Mal_vs_Ben = Malignant-Benign, levels=design)
fit2 <- eBayes(contrasts.fit(fit, cm))
tt <- topTable(fit2, coef=1, number=Inf, p.value=1)
tt$probe <- rownames(tt); tt$gene <- anno[tt$probe,"UCSC_RefGene_Name"]
write.csv(tt, file.path(TAB, "gse179458_DMP_Mal_vs_Ben.csv"), row.names=FALSE)
cat(sprintf("   DMP Malignant vs Benign: %d 显著 (FDR<0.05)\n", sum(tt$adj.P.Val<0.05)))

cat(">> DONE: 03_meth_gse179458.R\n")
