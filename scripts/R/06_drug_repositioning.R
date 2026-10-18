#!/usr/bin/env Rscript
# ============================================================
# 06_drug_repositioning.R — 候选基因 -> 可药靶网络 (PT 机制相关)
# 输入: 02/04/05 产出的显著基因 (DMP/DMR/DE)
# 方法: 收候选基因 -> DGIdb 可药靶注释(网络) + 回退内置关键通路靶标 -> PT机制分类 -> 乳腺表达过滤
# 输出: results/tables/gse_pt_druggable_candidates.csv
# ============================================================
suppressMessages({ library(limma); library(org.Hs.eg.db); library(httr); library(jsonlite); library(GEOquery) })

# 项目路径唯一真相源 (对齐 scripts/_env.py)
.cmdArgs <- commandArgs(trailingOnly = FALSE)
.fileArg <- sub("^--file=", "", .cmdArgs[grepl("^--file=", .cmdArgs)])
.sdir <- if (length(.fileArg)) normalizePath(dirname(.fileArg[1]), winslash="/") else normalizePath(getwd(), winslash="/")
source(file.path(.sdir, "_env.R"))

DS  <- DATASET_PATH
RAW <- file.path(DS, "raw")
PRC <- file.path(DS, "data/processed"); dir.create(PRC, recursive=TRUE)
RES <- file.path(DS, "results"); TAB <- file.path(RES, "tables"); dir.create(TAB, recursive=TRUE)
CACHE <- file.path(PRC, "dgidb_interactions.rds")

cat(">> [1] 收集候选基因 (02 DMP/DMR + 04 DE + 05 DE)\n")
cand <- c()
# 02 DMP (Mal_vs_Ben)
f <- file.path(TAB, "gse231574_DMP_Mal_vs_Ben.csv")
if (file.exists(f)) {
  d <- read.csv(f); sig <- d[d$adj.P.Val<0.05, "gene"]
  cand <- c(cand, unlist(strsplit(as.character(sig), ";")))
}
# 02 DMR (提取区域基因, 若 DMR 表含 gene 列)
f <- file.path(TAB, "gse231574_DMR_Mal_vs_Ben.csv")
if (file.exists(f)) {
  d <- read.csv(f)
  if ("gene" %in% colnames(d)) cand <- c(cand, unlist(strsplit(as.character(d$gene), ";")))
}
# 04 DE (Malignant vs Benign)
f <- file.path(TAB, "gse235267_DE_Mal_vs_Ben.csv")
if (file.exists(f)) { d <- read.csv(f); cand <- c(cand, as.character(d[d$adj.P.Val<0.05, "probe"])) }
# 05 DE (PT_vs_FA, Mal_vs_Ben) — 用 symbol 列
for (cn in c("PT_vs_FA","Mal_vs_Ben")) {
  f <- file.path(TAB, sprintf("gse78071_DE_%s.csv", cn))
  if (file.exists(f)) { d <- read.csv(f); s <- na.omit(d[d$adj.P.Val<0.05, "symbol"]); cand <- c(cand, as.character(s)) }
}
cand <- unique(cand[cand!="" & !is.na(cand)])
cat(sprintf("   候选基因合计 %d 个\n", length(cand)))

cat(">> [2] PT 机制相关通路靶标 (用于分类/回退)\n")
pt_pathways <- list(
  CellCycle_CDK = c("CDK1","CDK2","CDK4","CDK6","CCND1","CCNE1","RB1","CDKN2A"),
  PI3K_mTOR    = c("PIK3CA","PIK3CB","AKT1","AKT2","MTOR","PTEN","PDPK1"),
  RTK_PDGFR_FGFR = c("PDGFA","PDGFB","PDGFRA","PDGFRB","FGFR1","FGFR2","FGFR3","KIT","KDR","EGFR","ERBB2"),
  TGFb_ECM     = c("TGFB1","TGFB2","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","MMP2","MMP9","MMP14","COL1A1","CTGF"),
  Stromal_other= c("PDGFD","LOX","LOXL2","SNAI1","TWIST1","ZEB1","VIM","ACTA2","MCAM")
)
pt_genes <- unique(unlist(pt_pathways))

cat(">> [3] DGIdb 可药靶注释 (网络优先, 失败回退内置清单)\n")
drug_hits <- list()
if (!file.exists(CACHE)) {
  tryCatch({
    url <- paste0("https://dgidb.org/api/v1/interactions.json?genes=", paste(cand, collapse=","))
    resp <- GET(url, timeout(120))
    if (status_code(resp)==200) {
      j <- fromJSON(content(resp, "text", encoding="UTF-8"))
      # DGIdb 返回 list with $data$data 或 $interactions
      dat <- if (!is.null(j$data)) j$data else j
      saveRDS(dat, CACHE)
      drug_hits$dgidb <- dat
      cat("   DGIdb API OK, 缓存\n")
    } else stop("http ", status_code(resp))
  }, error = function(e) cat("   DGIdb API 失败:", conditionMessage(e), "-> 回退内置清单\n"))
} else {
  drug_hits$dgidb <- readRDS(CACHE); cat("   DGIdb 缓存已读\n")
}

# 解析 DGIdb: 每基因的药物数 + 类别
gene2drugs <- list()
if (!is.null(drug_hits$dgidb)) {
  d <- drug_hits$dgidb
  # DGIdb interactions.json 结构: $interactions 为 data.frame, 含 gene_name / drug_name / interaction_types / categories
  if (!is.null(d$interactions)) {
    it <- d$interactions
    for (g in unique(it$gene_name)) {
      sub <- it[it$gene_name==g, ]
      gene2drugs[[g]] <- list(n_drugs=nrow(sub),
                               drugs=paste(unique(sub$drug_name), collapse=";"),
                               cats=paste(unique(unlist(sub$categories)), collapse=";"))
    }
  }
}

cat(">> [4] 合并: 可药靶判定 + PT机制标签 + 乳腺表达过滤\n")
# 回退: 若 DGIdb 无命中, 用内置 pt_genes 判定"机制相关可药靶"
rows <- list()
for (g in cand) {
  dgidb_hit <- !is.null(gene2drugs[[g]])
  n_drug <- if (dgidb_hit) gene2drugs[[g]]$n_drugs else 0
  drugs  <- if (dgidb_hit) gene2drugs[[g]]$drugs else ""
  cats   <- if (dgidb_hit) gene2drugs[[g]]$cats else ""
  # PT 机制标签
  pw <- names(pt_pathways)[sapply(pt_pathways, function(v) g %in% v)]
  pt_tag <- if (length(pw)>0) paste(pw, collapse=";") else (if (g %in% pt_genes) "PT_other" else "")
  druggable <- dgidb_hit || (g %in% pt_genes)
  rows[[g]] <- data.frame(gene=g, druggable=druggable, n_drugs=n_drug, drugs=drugs,
                          dgidb_categories=cats, pt_pathway=pt_tag, stringsAsFactors=FALSE)
}
if (length(rows) == 0) {
  cand_df <- data.frame(gene=character(), druggable=logical(), n_drugs=integer(),
                        drugs=character(), dgidb_categories=character(),
                        pt_pathway=character(), stringsAsFactors=FALSE)
} else {
  cand_df <- do.call(rbind, rows)
}

# 乳腺表达过滤: 用 GSE78071 (Illumina HT-12) 表达均值 (PT vs FA)
fe <- file.path(RAW, "GSE78071/GSE78071_series_matrix.txt.gz")
if (file.exists(fe)) {
  lines <- readLines(gzfile(fe), warn=FALSE)
  beg <- grep("^!series_matrix_table_begin", lines); end <- grep("^!series_matrix_table_end", lines)
  m <- as.matrix(read.delim(text=lines[(beg+1):(end-1)], header=TRUE, check.names=FALSE, row.names=1))
  gpl <- tryCatch({ if (!file.exists(file.path(PRC,"GPL14951.rds"))) { g<-getGEO("GPL14951",destdir=PRC,AnnotGPL=TRUE); saveRDS(g,file.path(PRC,"GPL14951.rds")) } else readRDS(file.path(PRC,"GPL14951.rds")) }, error=function(e) NULL)
  if (!is.null(gpl)) {
    tab <- Table(gpl)
    idc <- if ("ID" %in% colnames(tab)) "ID" else colnames(tab)[1]
    symc <- grep("symbol", colnames(tab), ignore.case=TRUE)[1]
    if (is.na(symc)) symc <- grep("gene", colnames(tab), ignore.case=TRUE)[1]
    sm <- setNames(as.character(tab[[symc]]), as.character(tab[[idc]]))
    # 表型
    cl <- lines[grep("^!Sample_characteristics_ch1", lines)]
    cll <- lapply(cl, function(l) gsub('"','',unlist(strsplit(l,"\t")))[-1])
    nm <- max(sapply(cll,length)); cm <- do.call(cbind, lapply(cll, function(x) c(x,rep(NA,nm-length(x)))))
    sc <- apply(cm,1,function(x) paste(x[!is.na(x)&x!=""],collapse=" | "))
    grp <- ifelse(grepl("fibroadenoma",sc,ignore.case=TRUE),"FA",
          ifelse(grepl("benign phyllodes",sc,ignore.case=TRUE),"Benign_PT",
          ifelse(grepl("borderline phyllodes",sc,ignore.case=TRUE),"Borderline_PT",
          ifelse(grepl("malignant phyllodes",sc,ignore.case=TRUE),"Malignant_PT","Other"))))
    names(grp) <- colnames(m)
    for (g in cand_df$gene) {
      pr <- names(sm)[sm==g]
      if (length(pr)>0) {
        ex <- as.numeric(m[pr[1], ])
        cand_df[cand_df$gene==g, "mean_expr_PT"] <- mean(ex[grp %in% c("Benign_PT","Borderline_PT","Malignant_PT")], na.rm=TRUE)
        cand_df[cand_df$gene==g, "mean_expr_FA"] <- mean(ex[grp=="FA"], na.rm=TRUE)
      }
    }
    cat("   乳腺表达过滤完成 (GSE78071)\n")
  }
}

# 排序: 优先 PT机制相关 + 可药靶 + 药物多
cand_df$pt_rel <- cand_df$pt_pathway!="" & !is.na(cand_df$pt_pathway)
cand_df <- cand_df[order(-cand_df$pt_rel, -cand_df$druggable, -cand_df$n_drugs), ]
out <- file.path(TAB, "gse_pt_druggable_candidates.csv")
write.csv(cand_df, out, row.names=FALSE)
cat(sprintf(">> DONE: 06 可药靶候选 %d 个 (其中PT机制相关 %d, 有DGIdb药物 %d)\n   输出: %s\n",
            nrow(cand_df), sum(cand_df$pt_rel), sum(cand_df$n_drugs>0), out))
