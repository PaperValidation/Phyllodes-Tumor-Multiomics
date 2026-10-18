#!/usr/bin/env Rscript
suppressMessages(library(dplyr))
# 项目路径唯一真相源 (对齐 scripts/_env.py)
.cmdArgs <- commandArgs(trailingOnly = FALSE)
.fileArg <- sub("^--file=", "", .cmdArgs[grepl("^--file=", .cmdArgs)])
.sdir <- if (length(.fileArg)) normalizePath(dirname(.fileArg[1]), winslash="/") else normalizePath(getwd(), winslash="/")
source(file.path(.sdir, "_env.R"))

TAB <- TABLES
in_csv <- file.path(TAB, "gse_pt_druggable_candidates.csv")
out_csv <- file.path(TAB, "gse_pt_druggable_candidates.csv")
df <- read.csv(in_csv, stringsAsFactors=FALSE)
cat(sprintf("read %d candidates\n", nrow(df)))
pairs <- c(
 "CDK4|palbociclib;ribociclib;abemaciclib|CellCycle_CDK",
 "CDK6|palbociclib;ribociclib;abemaciclib|CellCycle_CDK",
 "CDK2|dinaciclib;seliciclib|CellCycle_CDK",
 "CDK1|dinaciclib;roscovitine|CellCycle_CDK",
 "CCND1|palbociclib(combo);ribociclib(combo)|CellCycle_CDK",
 "CCNE1|CDK4/6 inhibitors(investigational)|CellCycle_CDK",
 "CDKN2A|CDK4/6 inhibitors(context)|CellCycle_CDK",
 "RB1|CDK4/6 inhibitors(synthetic-lethal)|CellCycle_CDK",
 "PIK3CA|alpelisib;taselisib;copanlisib|PI3K_mTOR",
 "AKT1|capivasertib;ipatasertib|PI3K_mTOR",
 "MTOR|everolimus;temsirolimus;sirolimus|PI3K_mTOR",
 "PTEN|everolimus(context)|PI3K_mTOR",
 "PDGFRA|imatinib;avapritinib;sunitinib|RTK_PDGFR_FGFR",
 "PDGFRB|imatinib;sunitinib;pazopanib|RTK_PDGFR_FGFR",
 "PDGFA|imatinib;nintedanib|RTK_PDGFR_FGFR",
 "PDGFB|imatinib;nintedanib|RTK_PDGFR_FGFR",
 "PDGFD|PDGFR inhibitors(context)|RTK_PDGFR_FGFR",
 "FGFR1|erdafitinib;pemigatinib;infigratinib|RTK_PDGFR_FGFR",
 "FGFR2|erdafitinib;pemigatinib;infigratinib|RTK_PDGFR_FGFR",
 "FGFR3|erdafitinib;pemigatinib|RTK_PDGFR_FGFR",
 "KIT|imatinib;sunitinib;avapritinib|RTK_PDGFR_FGFR",
 "KDR|sunitinib;pazopanib;regorafenib;axitinib|RTK_PDGFR_FGFR",
 "EGFR|gefitinib;erlotinib;cetuximab;lapatinib|RTK_PDGFR_FGFR",
 "ERBB2|trastuzumab;lapatinib;pertuzumab;ado-trastuzumab|RTK_PDGFR_FGFR",
 "TGFB1|fresolimumab;galunisertib(investigational)|TGFb_ECM",
 "TGFB2|galunisertib(investigational)|TGFb_ECM",
 "TGFBR1|galunisertib;vactosertib(investigational)|TGFb_ECM",
 "TGFBR2|galunisertib(investigational)|TGFb_ECM",
 "MMP2|marimastat;prinomastat(investigational)|TGFb_ECM",
 "MMP9|doxycycline;marimastat(investigational)|TGFb_ECM",
 "MMP14|MMRI(investigational)|TGFb_ECM",
 "COL1A1|losartan(antifibrotic)|TGFb_ECM",
 "CTGF|pamrevlumab;FG-3019(investigational)|TGFb_ECM",
 "LOX|beta-aminopropionitrile;LOXL2 inhibitors(investigational)|Stromal_other",
 "LOXL2|simtuzumab(investigational);PXS-LOX2|Stromal_other",
 "SNAI1|no approved(targeted)|Stromal_other",
 "TWIST1|investigational(EMT)|Stromal_other",
 "ZEB1|no approved(targeted)|Stromal_other",
 "VIM|no approved|Stromal_other",
 "ACTA2|pirfenidone(antifibrotic)|Stromal_other",
 "MCAM|investigational|Stromal_other",
 "TP53|eprenetapopt(APR-246, investigational)|Tumor_suppressor",
 "MDM2|idasanutlin;nutlin-3(investigational)|Tumor_suppressor",
 "PARP1|olaparib;niraparib;talazoparib|DNA_repair",
 "BRCA1|olaparib;niraparib|DNA_repair",
 "BRCA2|olaparib;niraparib|DNA_repair",
 "ESR1|fulvestrant;elacestrant|Nuclear_receptor",
 "AR|enzalutamide;abiraterone;darolutamide|Nuclear_receptor",
 "MET|crizotinib;capmatinib;tepotinib|RTK_OTHER",
 "ALK|crizotinib;alectinib;lorlatinib|RTK_OTHER",
 "BRAF|vemurafenib;dabrafenib;encorafenib|MAPK",
 "MEK1|trametinib;binimetinib;selumetinib|MAPK",
 "MEK2|trametinib;binimetinib|MAPK"
)
dict_list <- lapply(strsplit(pairs, "|", fixed=TRUE), function(v) list(gene=v[1], drugs=v[2], cat=v[3]))
names(dict_list) <- sapply(dict_list, function(x) x$gene)
for (i in seq_len(nrow(df))) {
  g <- df$gene[i]
  if (g %in% names(dict_list)) {
    d <- dict_list[[g]]
    df$drugs[i] <- ifelse(is.na(df$drugs[i]) || df$drugs[i]=="", d$drugs, paste(unique(c(strsplit(df$drugs[i],";", fixed=TRUE)[[1]], strsplit(d$drugs,";", fixed=TRUE)[[1]])), collapse=";"))
    df$dgidb_categories[i] <- ifelse(is.na(df$dgidb_categories[i]) || df$dgidb_categories[i]=="", d$cat, paste(unique(c(strsplit(df$dgidb_categories[i],";", fixed=TRUE)[[1]], d$cat)), collapse=";"))
    df$druggable[i] <- TRUE
  }
}
df$n_drugs <- sapply(strsplit(as.character(df$drugs), ";", fixed=TRUE), function(x) sum(x!="" & !is.na(x)))
df$pt_rel <- df$pt_pathway!="" & !is.na(df$pt_pathway)
df <- df[order(-df$pt_rel, -df$druggable, -df$n_drugs), ]
write.csv(df, out_csv, row.names=FALSE)
nd <- sum(df$n_drugs>0)
cat(sprintf("DONE: annotated %d candidates; %d have >=1 curated drug\n", nrow(df), nd))
