suppressMessages(library(GEOquery))
# 项目路径唯一真相源 (对齐 scripts/_env.py)
.cmdArgs <- commandArgs(trailingOnly = FALSE)
.fileArg <- sub("^--file=", "", .cmdArgs[grepl("^--file=", .cmdArgs)])
.sdir <- if (length(.fileArg)) normalizePath(dirname(.fileArg[1]), winslash="/") else normalizePath(getwd(), winslash="/")
source(file.path(.sdir, "_env.R"))

PRC <- PROCESSED
out <- file.path(PRC, "GPL14951.rds")
softp <- file.path(PRC, "GPL14951.soft.gz")
softp2 <- file.path(PRC, "GPL14951.soft")
ok <- FALSE
for (i in 1:6) {
  unlink(c(softp, softp2))   # force fresh download (avoid reusing partial)
  gpl <- tryCatch(getGEO("GPL14951", destdir=PRC, AnnotGPL=TRUE),
                  error=function(e) { cat("attempt", i, "failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(gpl)) {
    saveRDS(gpl, out); ok <- TRUE
    cat("OK saved GPL14951.rds (attempt", i, ")\n"); break
  }
  Sys.sleep(5)
}
if (!ok) cat("ALL ATTEMPTS FAILED\n")
