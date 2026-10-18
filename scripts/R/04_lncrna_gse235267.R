#!/usr/bin/env Rscript
# ============================================================
# 04_lncrna_gse235267.R — lncRNA-mRNA 网络 (GSE235267, 16例)
# 表型: fibroadenoma(4) / benign PT(6) / malignant PT(6)  [无borderline]
# lncRNA 识别: GPL33505 平台 DESCRIPTION 字段关键词
# 注释: getGEO + Table (GEOquery 2.70 已移除 getGPL/getGSE; GPL=S4, 经 Table() 取注释表)
# ============================================================
suppressMessages({ library(limma); library(org.Hs.eg.db); library(GEOquery) })
# 项目路径唯一真相源 (对齐 scripts/_env.py)
.cmdArgs <- commandArgs(trailingOnly = FALSE)
.fileArg <- sub("^--file=", "", .cmdArgs[grepl("^--file=", .cmdArgs)])
.sdir <- if (length(.fileArg)) normalizePath(dirname(.fileArg[1]), winslash="/") else normalizePath(getwd(), winslash="/")
source(file.path(.sdir, "_env.R"))
source(file.path(.sdir, "enrich_helpers.R"))

DS  <- DATASET_PATH
RAW <- file.path(DS, "raw/GSE235267")
PRC <- file.path(DS, "data/processed"); dir.create(PRC, recursive=TRUE)
RES <- file.path(DS, "results"); TAB <- file.path(RES, "tables"); dir.create(TAB, recursive=TRUE)
GPL_CACHE <- file.path(PRC, "GPL33505.rds")

cat(">> [1] 解析 matrix (数据段)\n")
lines <- readLines(gzfile(file.path(RAW, "GSE235267_series_matrix.txt.gz")), warn=FALSE)
beg <- grep("^!series_matrix_table_begin", lines); end <- grep("^!series_matrix_table_end", lines)
dat <- as.matrix(read.delim(text=lines[(beg+1):(end-1)], header=TRUE, check.names=FALSE, row.names=1))

cat(">> [2] 表型 (稳健: 合并每样本全部 characteristics_ch1)\n")
char_lines <- lines[grep("^!Sample_characteristics_ch1", lines)]
char_list <- lapply(char_lines, function(l) gsub('"','',unlist(strsplit(l,"\t")))[-1])
ncol_max <- max(sapply(char_list, length))
char_mat <- do.call(cbind, lapply(char_list, function(x) c(x, rep(NA, ncol_max-length(x)))))
sample_chars <- apply(char_mat, 1, function(x) paste(x[!is.na(x)&x!=""], collapse=" | "))
names(sample_chars) <- colnames(dat); sample_chars <- sample_chars[colnames(dat)]
group <- ifelse(grepl("fibroadenoma", sample_chars, ignore.case=TRUE), "FA",
         ifelse(grepl("benign", sample_chars, ignore.case=TRUE) & grepl("phyllodes", sample_chars, ignore.case=TRUE), "Benign_PT",
         ifelse(grepl("malignant", sample_chars, ignore.case=TRUE) & grepl("phyllodes", sample_chars, ignore.case=TRUE), "Malignant_PT",
         ifelse(grepl("phyllodes", sample_chars, ignore.case=TRUE), "Borderline_PT", "Other"))))
cat("   分组: "); print(table(group))
writeLines(sample_chars, file.path(PRC, "gse235267_sample_chars.txt"))

cat(">> [3] GPL33505 注释 -> lncRNA (DESCRIPTION 关键词)\n")
if (!file.exists(GPL_CACHE)) {
  gpl <- getGEO("GPL33505", destdir=PRC, AnnotGPL=TRUE); saveRDS(gpl, GPL_CACHE)
} else {
  gpl <- readRDS(GPL_CACHE)
}
tab <- Table(gpl)
desc <- setNames(as.character(tab$DESCRIPTION), as.character(tab$ID))
ncrna_pat <- "non-protein coding RNA|antisense RNA|long intergenic|microRNA|snRNA|snoRNA|lincRNA|lncRNA|ribozyme|sRNA|RNA, (?:small|long)"
is_nc  <- function(d) grepl(ncrna_pat, d, ignore.case=TRUE)
is_mrna<- function(d) grepl("mRNA", d, ignore.case=TRUE) & !is_nc(d)
lnc_idx  <- which(is_nc(desc[rownames(dat)]))
mrna_idx <- which(is_mrna(desc[rownames(dat)]))
cat(sprintf("   识别 lncRNA 探针 %d / mRNA 探针 %d (共 %d)\n", length(lnc_idx), length(mrna_idx), nrow(dat)))

cat(">> [4] 差异表达 (Malignant_PT vs Benign_PT, 及 PT vs FA)\n")
design <- model.matrix(~0+factor(group, levels=c("FA","Benign_PT","Malignant_PT")))
colnames(design) <- c("FA","Benign_PT","Malignant_PT")
fit <- lmFit(dat, design)
cm <- makeContrasts(Mal_vs_Ben = Malignant_PT-Benign_PT,
                    PT_vs_FA   = (Benign_PT+Malignant_PT)/2 - FA, levels=design)
fit2 <- eBayes(contrasts.fit(fit, cm))
for (cn in colnames(cm)) {
  tt <- topTable(fit2, coef=cn, number=Inf, p.value=1)
  tt$probe <- rownames(tt)
  tt$is_lncRNA <- rownames(tt) %in% rownames(dat)[lnc_idx]
  write.csv(tt, file.path(TAB, sprintf("gse235267_DE_%s.csv", cn)), row.names=FALSE)
  cat(sprintf("   DE %-12s: %d 显著 | 其中lncRNA %d\n", cn,
              sum(tt$adj.P.Val<0.05), sum(tt$adj.P.Val<0.05 & tt$is_lncRNA)))
}

cat(">> [5] lncRNA-mRNA 相关网络 (恶性/PT 相关 DE lncRNA, Spearman Top50)\n")
lnc_de <- c()
for (cn in c("Mal_vs_Ben","PT_vs_FA")) {
  df <- read.csv(file.path(TAB, sprintf("gse235267_DE_%s.csv", cn)))
  lnc_de <- unique(c(lnc_de, as.character(df[df$adj.P.Val<0.05 & df$is_lncRNA, "probe"])))
}
if (length(lnc_de) < 5) {   # 小队列 FDR 显著 lncRNA 极少 -> 放宽取 PT_vs_FA 中 |logFC| 最大的 lncRNA
  df <- read.csv(file.path(TAB, "gse235267_DE_PT_vs_FA.csv"))
  ln <- df[df$is_lncRNA, ]
  if (nrow(ln) > 0) lnc_de <- unique(c(lnc_de, as.character(ln$probe[order(-abs(ln$logFC))])[1:20]))
}
cat(sprintf("   网络种子 lncRNA %d 个\n", length(lnc_de)))
net <- NULL
if (length(lnc_de) > 0) {
  edges <- list()
  for (lp in lnc_de) {
    Xt <- t(dat[lp, , drop=FALSE])                       # n_samples x 1 (lncRNA)
    Yt <- t(dat[mrna_idx, , drop=FALSE])                 # n_samples x m (mRNA)
    r <- cor(Xt, Yt, method="spearman")                 # 1 x m correlations
    rvec <- as.numeric(r[1, ]); names(rvec) <- colnames(r)
    ord <- order(-abs(rvec))[1:50]
    edges[[lp]] <- data.frame(lncRNA=lp, mRNA=names(rvec)[ord], rho=as.numeric(rvec[ord]), stringsAsFactors=FALSE)
  }
  net <- do.call(rbind, edges)
  write.csv(net, file.path(TAB, "gse235267_lncRNA_mRNA_network.csv"), row.names=FALSE)
  cat(sprintf("   网络边 %d (涉及 %d 个 lncRNA)\n", nrow(net), length(lnc_de)))
} else {
  write.csv(data.frame(lncRNA=character(), mRNA=character(), rho=numeric()),
            file.path(TAB, "gse235267_lncRNA_mRNA_network.csv"), row.names=FALSE)
  cat("   (无可用 lncRNA 种子, 网络为空)\n")
}

cat(">> [6] mRNA伙伴功能富集 (stromal proliferation / cell cycle / ECM)\n")
extract_sym <- function(d) {
  m <- regmatches(d, regexpr("\\(([A-Za-z0-9_-]+)\\)", d))
  if (length(m) && !is.na(m) && m!="") return(gsub("[()]","",m))
  m2 <- regmatches(d, regexpr("Homo sapiens ([A-Za-z0-9_-]+)", d))
  if (length(m2) && !is.na(m2) && m2!="") return(gsub("Homo sapiens ","",m2))
  NA
}
if (is.null(net) || nrow(net)==0) {
  write.csv(data.frame(term=character(), pvalue=numeric(), stringsAsFactors=FALSE),
            file.path(TAB, "gse235267_lncRNA_partner_GO.csv"), row.names=FALSE)
  cat("   (网络为空, 跳过 GO 富集)\n")
} else {
  mrna_partners <- unique(net$mRNA)
  sym <- sapply(desc[mrna_partners], extract_sym)
  sym <- sym[!is.na(sym) & sym!=""]
  eg <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys=sym, columns="ENTREZID", keytype="SYMBOL"),
                 error=function(e) data.frame(SYMBOL=character(), ENTREZID=character()))
  eg <- eg[!is.na(eg$ENTREZID) & eg$ENTREZID!="", ]
  if (nrow(eg)==0) {
    write.csv(data.frame(term=character(), pvalue=numeric(), stringsAsFactors=FALSE),
              file.path(TAB, "gse235267_lncRNA_partner_GO.csv"), row.names=FALSE)
    cat(sprintf("   (mRNA伙伴 %d, 但无符号映射, 跳过 GO)\n", length(mrna_partners)))
  } else {
    go <- enrichGO_manual(eg$ENTREZID, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.1)
    go.df <- as.data.frame(go)
    write.csv(go.df, file.path(TAB, "gse235267_lncRNA_partner_GO.csv"), row.names=FALSE)
    cat(sprintf("   GO 条目 %d (输入mRNA伙伴 %d, 映射符号 %d)\n", nrow(go.df), length(mrna_partners), length(sym)))
  }
}

cat(">> DONE: 04_lncrna_gse235267.R\n")
