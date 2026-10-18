#!/usr/bin/env Rscript
# ============================================================
# enrich_helpers.R — 零依赖 GO 富集 (替代 clusterProfiler::enrichGO)
# 算法: 超几何检验 phyper, 背景 = org.Hs.eg.db 中全部带 GO 注释的 ENTREZ 基因
# 用法: source("enrich_helpers.R"); go <- enrichGO_manual(entrez_ids, ont="BP")
# 输出: data.frame(GO, Term, Count, GeneRatio, BgRatio, pvalue, p.adjust)
# ============================================================
suppressMessages({ library(org.Hs.eg.db); library(GO.db); library(AnnotationDbi) })

enrichGO_manual <- function(gene, universe=NULL, ont="BP",
                            pAdjustMethod="BH", pvalueCutoff=0.05,
                            minGSSize=10, maxGSSize=500, qvalueCutoff=0.2) {
  gene <- unique(as.character(gene[!is.na(gene) & gene != ""]))
  egGO <- org.Hs.egGO
  if (is.null(universe)) universe <- ls(egGO)
  universe <- intersect(unique(as.character(universe)), ls(egGO))
  N <- length(universe); n <- length(gene)
  if (n == 0 || N == 0) return(data.frame())

  ont_want <- if (toupper(ont) == "ALL") c("BP", "CC", "MF") else strsplit(ont, "/")[[1]]

  # 标准取 gene -> GO -> ontology (AnnotationDbi::select)
  anno <- tryCatch({
    AnnotationDbi::select(org.Hs.eg.db, keys = universe,
                          columns = c("GO", "ONTOLOGY"), keytype = "ENTREZID")
  }, error = function(e) data.frame(ENTREZID = character(), GO = character(), ONTOLOGY = character()))
  if (nrow(anno) == 0) return(data.frame())
  anno <- anno[!is.na(anno$GO) & anno$ONTOLOGY %in% ont_want, c("ENTREZID", "GO")]
  anno <- anno[!duplicated(anno), ]
  term_genes <- split(anno$ENTREZID, anno$GO)

  gene_set <- gene
  rows <- list()
  for (t in names(term_genes)) {
    K <- length(unique(term_genes[[t]]))
    if (K < minGSSize || K > maxGSSize) next
    k <- length(intersect(gene_set, unique(term_genes[[t]])))
    if (k == 0) next
    p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    rows[[t]] <- data.frame(GO = t, Count = k, pvalue = p, stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) return(data.frame())
  res <- do.call(rbind, rows)
  res$p.adjust <- p.adjust(res$pvalue, method = pAdjustMethod)
  res <- res[res$p.adjust <= qvalueCutoff, , drop = FALSE]
  if (nrow(res) == 0) return(data.frame())

  go_ids <- as.character(res$GO)
  desc <- AnnotationDbi::select(GO.db, keys = go_ids, columns = "TERM", keytype = "GOID")
  desc <- desc[!duplicated(desc$GOID), ]
  res <- merge(res, desc, by.x = "GO", by.y = "GOID", all.x = TRUE)
  if ("TERM" %in% colnames(res)) colnames(res)[colnames(res) == "TERM"] <- "Term"
  res$GeneRatio <- sprintf("%d/%d", res$Count, n)
  res$BgRatio   <- sapply(res$GO, function(t) length(unique(term_genes[[t]])))
  res$BgRatio   <- sprintf("%d/%d", res$BgRatio, N)
  res <- res[order(res$p.adjust), ]; rownames(res) <- NULL
  res
}
