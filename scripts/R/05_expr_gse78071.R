#!/usr/bin/env Rscript
# ============================================================
# 05_expr_gse78071.R — 表达验证 (GSE78071, 48例 FA/PT 芯片)
# 表型(pathological diagnosis): Simple/Complex fibroadenoma,
#   Benign/Borderline/Malignant phyllodes tumor  [完整进展谱系]
# 注释: GPL14951 (Illumina HumanHT-12) via getGEO + Table -> Gene.Symbol
# 流程: 加载matrix -> 表型 -> DE -> 甲基化模块 concordance(符号级) -> 恶性样程序异质性
# ============================================================
suppressMessages({ library(limma); library(org.Hs.eg.db); library(GEOquery) })

# 项目路径唯一真相源 (对齐 scripts/_env.py)
.cmdArgs <- commandArgs(trailingOnly = FALSE)
.fileArg <- sub("^--file=", "", .cmdArgs[grepl("^--file=", .cmdArgs)])
.sdir <- if (length(.fileArg)) normalizePath(dirname(.fileArg[1]), winslash="/") else normalizePath(getwd(), winslash="/")
source(file.path(.sdir, "_env.R"))

DS  <- DATASET_PATH
RAW <- file.path(DS, "raw/GSE78071")
PRC <- file.path(DS, "data/processed"); dir.create(PRC, recursive=TRUE)
RES <- file.path(DS, "results"); TAB <- file.path(RES, "tables"); dir.create(TAB, recursive=TRUE)
GPL_CACHE <- file.path(PRC, "GPL14951.rds")

cat(">> [1] 解析 matrix + 稳健表型\n")
lines <- readLines(gzfile(file.path(RAW, "GSE78071_series_matrix.txt.gz")), warn=FALSE)
beg <- grep("^!series_matrix_table_begin", lines); end <- grep("^!series_matrix_table_end", lines)
dat <- as.matrix(read.delim(text=lines[(beg+1):(end-1)], header=TRUE, check.names=FALSE, row.names=1))
char_lines <- lines[grep("^!Sample_characteristics_ch1", lines)]
char_list <- lapply(char_lines, function(l) gsub('"','',unlist(strsplit(l,"\t")))[-1])
ncol_max <- max(sapply(char_list, length))
char_mat <- do.call(cbind, lapply(char_list, function(x) c(x, rep(NA, ncol_max-length(x)))))
sample_chars <- apply(char_mat, 1, function(x) paste(x[!is.na(x)&x!=""], collapse=" | "))
names(sample_chars) <- colnames(dat); sample_chars <- sample_chars[colnames(dat)]
class <- ifelse(grepl("fibroadenoma", sample_chars, ignore.case=TRUE), "FA",
        ifelse(grepl("benign phyllodes", sample_chars, ignore.case=TRUE), "Benign_PT",
        ifelse(grepl("borderline phyllodes", sample_chars, ignore.case=TRUE), "Borderline_PT",
        ifelse(grepl("malignant phyllodes", sample_chars, ignore.case=TRUE), "Malignant_PT", "Other"))))
names(class) <- colnames(dat)   # FIX: attach sample IDs so class[names(score)] resolves correctly
cat("   分组: "); print(table(class))
writeLines(sample_chars, file.path(PRC, "gse78071_sample_chars.txt"))

cat(">> [2] GPL14951 注释 -> Gene.Symbol\n")
if (!file.exists(GPL_CACHE)) {
  gpl <- getGEO("GPL14951", destdir=PRC, AnnotGPL=TRUE); saveRDS(gpl, GPL_CACHE)
} else {
  gpl <- readRDS(GPL_CACHE)
}
tab <- Table(gpl)
idc <- if ("ID" %in% colnames(tab)) "ID" else colnames(tab)[1]
symc <- grep("symbol", colnames(tab), ignore.case=TRUE)[1]
if (is.na(symc)) symc <- grep("gene", colnames(tab), ignore.case=TRUE)[1]
sym_map <- setNames(as.character(tab[[symc]]), as.character(tab[[idc]]))
sym_map[is.na(sym_map) | sym_map==""] <- NA

cat(">> [3] 差异表达 (FA vs PT, 进展趋势, 恶性vs良性)\n")
design <- model.matrix(~0+factor(class, levels=c("FA","Benign_PT","Borderline_PT","Malignant_PT")))
colnames(design) <- c("FA","Benign_PT","Borderline_PT","Malignant_PT")

# --- FIX: Illumina HumanHT-12 标准化 ---
# GEO series matrix 不包含阴性对照探针 (GPL14951 的 Probe_Type 全为 A/I/M/S,
# 矩阵内 0 个 control 行), 因此 limma::neqc() 的 normexp 背景校正无法执行。
# 改用 log2 变换 + 分位数标准化 (neqc 的 "q" 步), 这是无对照探针时的标准做法,
# 远优于在 raw 强度上直接跑 lmFit (旧版 logFC 达 ±10000)。
dat <- log2(dat)
dat <- limma::normalizeBetweenArrays(dat, method="quantile")
cat(sprintf("   标准化后 log2 表达范围: %.2f ~ %.2f\n", min(dat), max(dat)))

fit <- lmFit(dat, design)
cm <- makeContrasts(PT_vs_FA = (Benign_PT+Borderline_PT+Malignant_PT)/3 - FA,
                    Trend    = (Borderline_PT-Benign_PT)+(Malignant_PT-Borderline_PT),
                    Mal_vs_Ben = Malignant_PT-Benign_PT, levels=design)
fit2 <- eBayes(contrasts.fit(fit, cm))
for (cn in colnames(cm)) {
  tt <- topTable(fit2, coef=cn, number=Inf, p.value=1)
  tt$probe <- rownames(tt)
  tt$symbol <- sym_map[tt$probe]
  write.csv(tt, file.path(TAB, sprintf("gse78071_DE_%s.csv", cn)), row.names=FALSE)
  cat(sprintf("   DE %-12s: %d 显著 (FDR<0.05) | 有符号映射 %d\n",
              cn, sum(tt$adj.P.Val<0.05), sum(tt$adj.P.Val<0.05 & !is.na(tt$symbol))))
}

cat(">> [4] 甲基化高风脸模块基因 -> 表达层 concordance (符号级)\n")
meth_dmp <- NULL
for (f in c("gse231574_DMP_Mal_vs_Ben.csv","gse179458_DMP_Mal_vs_Ben.csv")) {
  p <- file.path(TAB, f); if (file.exists(p)) meth_dmp <- unique(c(meth_dmp, read.csv(p)$gene))
}
meth_dmp <- unique(unlist(strsplit(as.character(meth_dmp), ";")))
meth_dmp <- meth_dmp[meth_dmp!="" & !is.na(meth_dmp)]
de_pt <- read.csv(file.path(TAB, "gse78071_DE_PT_vs_FA.csv"))
de_sym <- unique(na.omit(de_pt$symbol))
overlap <- intersect(meth_dmp, de_sym)
cat(sprintf("   甲基化DMP基因 %d 个, 与 GSE78071 PTvsFA DE(符号) 重叠 %d 个 (%.1f%%)\n",
            length(meth_dmp), length(overlap), 100*length(overlap)/max(1,length(meth_dmp))))
if (length(overlap)>0) {
  sub <- de_pt[!is.na(de_pt$symbol) & de_pt$symbol %in% overlap, c("symbol","logFC","adj.P.Val")]
  write.csv(unique(sub), file.path(TAB, "gse78071_meth_module_concordance.csv"), row.names=FALSE)
}

cat(">> [5] 恶性样程序异质性 (benign/borderline PT 中高表达恶性signature者)\n")
mal_sig <- read.csv(file.path(TAB, "gse78071_DE_Mal_vs_Ben.csv"))
up <- mal_sig[mal_sig$adj.P.Val<0.05 & mal_sig$logFC>0, "probe"]
if (length(up) >= 10) {
  score <- apply(dat[up,,drop=FALSE], 2, function(x) mean(x, na.rm=TRUE))
  sc_df <- data.frame(sample=names(score), score=as.numeric(score),
                      class=class[names(score)], stringsAsFactors=FALSE)
  write.csv(sc_df, file.path(TAB, "gse78071_malignant_signature_score.csv"), row.names=FALSE)
  mal_med <- median(score[class[names(score)]=="Malignant_PT"])
  hot <- sc_df[sc_df$class %in% c("Benign_PT","Borderline_PT") & sc_df$score>mal_med, ]
  cat(sprintf("   恶性signature上调探针 %d; benign/borderline 中已呈恶性样程序样本 %d 例:\n", length(up), nrow(hot)))
  print(hot)
} else cat("   (恶性DE上调探针不足, 跳过异质性分析)\n")

cat(">> DONE: 05_expr_gse78071.R\n")
