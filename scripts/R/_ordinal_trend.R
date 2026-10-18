## =====================================================================
## _ordinal_trend.R  --  TRUE ordinal trend models (continuous severity)
## =====================================================================
## Motivation: the legacy "Trend" contrasts are algebraically degenerate:
##     (Borderline - Benign) + (Malignant - Borderline)  ===  Malignant - Benign
## so the legacy Trend tables are byte-identical to the pairwise tables and
## carried no independent evidence. This script refits a genuine ordinal
## (continuous-severity) model.
##
## Pipeline is VERIFIED against the legacy tables before use:
##   * GSE78071 PT_vs_FA  -> 737 probes, min P = 7.8959e-09  (legacy: 737, 7.8959e-09)
##   * GSE78071 Mal_vs_Ben-> 62 probes                          (legacy: 62)
##   * GSE231574 Mal_vs_Ben -> 0 probes, min adj.P = 1.000000, slope 1.000 (legacy identical)
##
## READ-ONLY w.r.t. all existing data/scripts. Writes ONLY new files:
##   results/tables/gse78071_DE_ordinal_trend.csv            (spec model, PT only)
##   results/tables/gse78071_DE_ordinal_trend_pooledvar.csv  (same slope, variance pooled over all 48)
##   results/tables/gse231574_DMP_ordinal_trend.csv
##   results/tables/_ordinal_trend_summary.txt
## =====================================================================
suppressPackageStartupMessages({
  library(limma)
  library(GEOquery)
})

## ---- project paths: single source of truth (.env) ----
.cmdArgs <- commandArgs(trailingOnly = FALSE)
.fileArg <- sub("^--file=", "", .cmdArgs[grepl("^--file=", .cmdArgs)])
.sdir <- if (length(.fileArg)) normalizePath(dirname(.fileArg[1]), winslash="/") else normalizePath(getwd(), winslash="/")
source(file.path(.sdir, "_env.R"))

ROOT <- DATASET_PATH
PRC  <- PROCESSED
TAB  <- TABLES
stopifnot(dir.exists(ROOT), dir.exists(PRC), dir.exists(TAB))

sink(file.path(TAB, "_ordinal_trend_summary.txt"), split = TRUE)
cat("=====================================================================\n")
cat(" TRUE ORDINAL TREND MODELS   ran:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(" severity: Benign=1 / Borderline=2 / Malignant=3 (continuous covariate)\n")
cat(" limma lmFit + eBayes + BH ; logFC = slope per grade step\n")
cat("=====================================================================\n\n")

## =====================================================================
## PART 1 -- GSE78071 expression (WG-DASL, GPL14951)
## =====================================================================
lines <- readLines(gzfile(file.path(ROOT, "raw/GSE78071/GSE78071_series_matrix.txt.gz")), warn = FALSE)
ch <- lines[startsWith(lines, "!Sample_characteristics_ch1")]
diag <- sub("^.*pathological diagnosis: *", "",
            gsub('^"|"$', "", strsplit(ch[grepl("pathological diagnosis", ch)][1], "\t")[[1]][-1]))
gsm <- gsub('^"|"$', "", strsplit(lines[startsWith(lines, "!Sample_geo_accession")][1], "\t")[[1]][-1])
tb <- grep("^!series_matrix_table_begin", lines); te <- grep("^!series_matrix_table_end", lines)
body <- lines[(tb + 1):(te - 1)]; hdr <- gsub('^"|"$', "", strsplit(body[1], "\t")[[1]])
raw <- as.matrix(data.table::fread(text = paste(body[-1], collapse = "\n"), header = FALSE,
                                   data.table = FALSE, showProgress = FALSE))
probe <- raw[, 1]; raw <- apply(raw[, -1, drop = FALSE], 2, as.numeric); rownames(raw) <- probe
colnames(raw) <- hdr[-1]
cat("GSE78071:", dim(raw), "samples | diagnosis:\n"); print(table(diag))

## normalisation: log2 then quantile.  log2(x) (NOT log2(x+1)) is what reproduces
## the legacy tables exactly (737, min P = 7.8959e-09); GPL14951 has no negative
## control probes, so neqc is not applicable.
E <- normalizeBetweenArrays(log2(raw), method = "quantile"); colnames(E) <- hdr[-1]

gpl <- readRDS(file.path(PRC, "GPL14951.rds")); anno <- Table(gpl)[, c("ID", "Symbol")]
sym <- anno$Symbol[match(rownames(E), anno$ID)]

cls <- ifelse(grepl("^Simple|^Complex", diag), "FA",
       ifelse(grepl("^Benign phyllodes", diag), "Ben",
       ifelse(grepl("^Borderline phyllodes", diag), "Bor", "Mal")))
grp <- factor(cls, levels = c("FA", "Ben", "Bor", "Mal"))
cat("\ngroup sizes:\n"); print(table(grp))
d4 <- model.matrix(~0 + grp); colnames(d4) <- levels(grp)

## ---- pipeline verification against legacy tables (not written) ----
cmV <- makeContrasts(PTvsFA = "(Ben+Bor+Mal)/3 - FA", levels = d4)
ttV <- topTable(eBayes(contrasts.fit(lmFit(E, d4), cmV)), coef = 1, number = Inf, sort.by = "P")
cat(sprintf("[verify] PT_vs_FA  = %d probes FDR<0.05 (legacy 737), min P = %.4e (legacy 7.8959e-09)\n",
            sum(ttV$adj.P.Val < 0.05), min(ttV$P.Value)))
cmV2 <- makeContrasts(MalvsBen = "Mal - Ben", levels = d4)
ttV2 <- topTable(eBayes(contrasts.fit(lmFit(E, d4), cmV2)), coef = 1, number = Inf, sort.by = "P")
cat(sprintf("[verify] Mal_vs_Ben = %d probes FDR<0.05 (legacy 62), min P = %.4e (legacy 2.6935e-10)\n\n",
            sum(ttV2$adj.P.Val < 0.05), min(ttV2$P.Value)))

## ---- MODEL A (spec): severity regression on the PT subset only ----
sub <- cls != "FA"
S <- as.integer(factor(cls[sub], levels = c("Ben", "Bor", "Mal")))
cat("[A] spec model: ~severity, PT subset only, n =", sum(sub), ":", paste(names(table(S)), table(S), collapse = " / "), "\n")
dA <- model.matrix(~severity, data = data.frame(severity = S))
fitA <- eBayes(lmFit(E[, sub], dA))
ttA <- topTable(fitA, coef = "severity", number = Inf, sort.by = "P")
outA <- data.frame(logFC = ttA$logFC, AveExpr = ttA$AveExpr, t = ttA$t, P.Value = ttA$P.Value,
                   adj.P.Val = ttA$adj.P.Val, B = ttA$B, probe = rownames(ttA),
                   symbol = sym[match(rownames(ttA), rownames(E))], stringsAsFactors = FALSE)
fA <- file.path(TAB, "gse78071_DE_ordinal_trend.csv")
write.csv(outA, fA, row.names = FALSE)
cat("WROTE:", fA, " rows:", nrow(outA), "\n")
cat(sprintf("  FDR<0.05 = %d probe(s) / %d gene(s) | nominal p<0.05 = %d (%.2fx expected) | p<0.001 = %d (expected %.0f)\n",
            sum(outA$adj.P.Val < 0.05),
            length(unique(outA$symbol[outA$adj.P.Val < 0.05])),
            sum(outA$P.Value < 0.05), sum(outA$P.Value < 0.05) / (0.05 * nrow(outA)),
            sum(outA$P.Value < 0.001), 0.001 * nrow(outA)))
cat("  df.total =", round(fitA$df.residual[1] + fitA$df.prior, 2), "| min adj.P =", signif(min(outA$adj.P.Val), 4), "\n")
print(utils::head(outA[order(outA$P.Value), c("probe", "symbol", "logFC", "t", "P.Value", "adj.P.Val")], 8))

## ---- MODEL B: identical slope, variance pooled over all 48 samples ----
## centred linear contrast over grades, FA contributes to the variance only:
## w_i = n_i (x_i - xbar) with x = (NA,1,2,3), n = (24,16,4,4); sum n_i(x_i-xbar)^2 = 14
## so contrast (0,-8,+2,+6)/14 returns the OLS slope PER GRADE STEP (same units as model A)
cmB <- makeContrasts(slope = (0*FA - 8*Ben + 2*Bor + 6*Mal)/14, levels = d4)
fitB <- eBayes(contrasts.fit(lmFit(E, d4), cmB))
ttB <- topTable(fitB, coef = 1, number = Inf, sort.by = "P")
outB <- data.frame(logFC = ttB$logFC, AveExpr = ttB$AveExpr, t = ttB$t, P.Value = ttB$P.Value,
                   adj.P.Val = ttB$adj.P.Val, B = ttB$B, probe = rownames(ttB),
                   symbol = sym[match(rownames(ttB), rownames(E))], stringsAsFactors = FALSE)
fB <- file.path(TAB, "gse78071_DE_ordinal_trend_pooledvar.csv")
write.csv(outB, fB, row.names = FALSE)
cat("\nWROTE:", fB, " rows:", nrow(outB), "\n")
cat(sprintf("[B] pooled variance (df.total = %.1f, same as legacy model): FDR<0.05 = %d probe(s) / %d gene(s) | nominal = %d (%.2fx) | p<0.001 = %d (%.1fx)\n",
            fitB$df.residual[1] + fitB$df.prior, sum(outB$adj.P.Val < 0.05),
            length(unique(outB$symbol[outB$adj.P.Val < 0.05])),
            sum(outB$P.Value < 0.05), sum(outB$P.Value < 0.05) / (0.05 * nrow(outB)),
            sum(outB$P.Value < 0.001), sum(outB$P.Value < 0.001) / (0.001 * nrow(outB))))
print(utils::head(outB[order(outB$P.Value), c("probe", "symbol", "logFC", "t", "P.Value", "adj.P.Val")], 15))

## ---- comparison with the legacy pairwise table ----
oldMB <- data.table::fread(file.path(TAB, "gse78071_DE_Mal_vs_Ben.csv"), data.table = FALSE, showProgress = FALSE)
cat("\n--- new ordinal models vs legacy Mal_vs_Ben (62 hits) ---\n")
for (tag in c("A", "B")) {
  tt <- if (tag == "A") ttA else ttB
  m <- merge(data.frame(probe = rownames(tt), lf = tt$logFC), oldMB[, c("probe", "logFC")], by = "probe")
  hits <- rownames(tt)[tt$adj.P.Val < 0.05]
  cat(sprintf("  [%s] Spearman(logFC) = %.4f | Pearson = %.4f | hits = %d | overlap with legacy 62 = %d | identical logFC vectors? %s\n",
              tag, cor(m$lf, m$logFC, method = "spearman"), cor(m$lf, m$logFC),
              length(hits), length(intersect(hits, oldMB$probe[oldMB$adj.P.Val < 0.05])),
              isTRUE(all.equal(m$lf, m$logFC))))
}
cat("\n  (legacy 'Trend' == Mal_vs_Ben by construction: (Bor-Ben)+(Mal-Bor) = Mal-Ben)\n")

## =====================================================================
## PART 2 -- GSE231574 methylation (M values)
## =====================================================================
cat("\n### PART 2: GSE231574 methylation ordinal trend ###\n")
M <- readRDS(file.path(PRC, "gse231574_Mval_filtered.rds"))
ph <- readRDS(file.path(PRC, "gse231574_pheno.rds"))
cat("M matrix:", dim(M), "| pheno:", dim(ph), "\n")
g5 <- factor(ph$group3, levels = c("FA", "Benign", "Borderline", "Malignant", "Metaplastic"))
names(g5) <- as.character(ph$Sample.ID)
print(table(g5))
col_g <- g5[match(colnames(M), names(g5))]
dM5 <- model.matrix(~0 + g5); colnames(dM5) <- levels(g5)

## verification: reproduce the legacy Mal_vs_Ben table (0 hits, min adj.P = 1)
cmM <- makeContrasts(MalvsBen = Malignant - Benign, levels = dM5)
ttV3 <- topTable(eBayes(contrasts.fit(lmFit(M, dM5), cmM)), coef = 1, number = Inf, sort.by = "P")
cat(sprintf("[verify] Mal_vs_Ben = %d probes FDR<0.05 (legacy 0), min adj.P = %.6f (legacy 0.99999957)\n",
            sum(ttV3$adj.P.Val < 0.05), min(ttV3$adj.P.Val)))

sevM <- as.integer(factor(as.character(col_g), levels = c("Benign", "Borderline", "Malignant")))
keepP <- !is.na(sevM)
cat("PT subset n =", sum(keepP), "(dropped", sum(!keepP), "FA/Metaplastic)\n")
print(table(as.character(col_g)[keepP]))
dM2 <- model.matrix(~severity, data = data.frame(severity = sevM[keepP]))
fitOM <- eBayes(lmFit(M[, keepP], dM2))
ttOM <- topTable(fitOM, coef = "severity", number = Inf, sort.by = "P")
ttOM$probe <- rownames(ttOM)
old <- data.table::fread(file.path(TAB, "gse231574_DMP_Mal_vs_Ben.csv"), data.table = FALSE, showProgress = FALSE)
ttOM$gene   <- old$gene[match(ttOM$probe, old$probe)]
ttOM$region <- old$region[match(ttOM$probe, old$probe)]
out2 <- ttOM[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "probe", "gene", "region")]
f2 <- file.path(TAB, "gse231574_DMP_ordinal_trend.csv")
write.csv(out2, f2, row.names = FALSE)
cat("WROTE:", f2, " rows:", nrow(out2), "\n")
cat(sprintf("[meth] FDR<0.05 = %d | nominal p<0.05 = %d (%.2fx expected) | p<0.001 = %d (%.2fx) | min adj.P = %.4f | df.total = %.1f\n",
            sum(out2$adj.P.Val < 0.05), sum(out2$P.Value < 0.05),
            sum(out2$P.Value < 0.05) / (0.05 * nrow(out2)),
            sum(out2$P.Value < 0.001), sum(out2$P.Value < 0.001) / (0.001 * nrow(out2)),
            min(out2$adj.P.Val), fitOM$df.residual[1] + fitOM$df.prior))
m2 <- merge(data.frame(probe = ttOM$probe, lf = ttOM$logFC),
            old[, c("probe", "logFC")], by = "probe")
cat(sprintf("  vs legacy Mal_vs_Ben: Spearman(logFC) = %.4f | hits_new = %d | overlap = %d | identical logFC? %s\n",
            cor(m2$lf, m2$logFC, method = "spearman"), sum(out2$adj.P.Val < 0.05),
            sum(out2$adj.P.Val < 0.05 & old$adj.P.Val[match(out2$probe, old$probe)] < 0.05),
            isTRUE(all.equal(m2$lf, m2$logFC))))

cat("\n### DONE ###\n")
sink()
