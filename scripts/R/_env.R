# ============================================================
# _env.R — R 侧项目路径唯一真相源（对齐 scripts/_env.py）
#
# 用途: 解析项目根 .env 的 DATASET_PATH / DATASET_PATH_TEMP,
#       并导出 RAW / PROCESSED / TABLES / FIGURES / TEMP 等路径变量,
#       使 R 脚本不再硬编码 /mnt/g/... 或 /mnt/c/...。
#
# 用法 (在每个 R 脚本顶部):
#   .cmdArgs <- commandArgs(trailingOnly = FALSE)
#   .fileArg <- sub("^--file=", "", .cmdArgs[grepl("^--file=", .cmdArgs)])
#   .sdir <- if (length(.fileArg)) normalizePath(dirname(.fileArg[1]), winslash="/")
#            else normalizePath(getwd(), winslash="/")
#   source(file.path(.sdir, "_env.R"))
#
# 之后即可使用 DATASET_PATH / RAW / PROCESSED / TABLES / FIGURES / TEMP。
# ============================================================

# ---------- 工具: 解析 .env 中单个键 ----------
.env_parse <- function(key) {
  lines <- tryCatch(
    readLines(file.path(.PT_PROJECT_ROOT, ".env"), warn = FALSE),
    error = function(e) character(0)
  )
  for (ln in lines) {
    ln <- sub("^[ \t]+", "", ln)
    if (!nzchar(ln) || startsWith(ln, "#") || !grepl("=", ln)) next
    kv <- strsplit(ln, "=", fixed = TRUE)[[1]]
    if (trimws(kv[1]) == key) {
      val <- trimws(paste(kv[-1], collapse = "="))
      # 去掉单/双引号
      val <- sub("^['\"]", "", val); val <- sub("['\"]$", "", val)
      return(trimws(val))
    }
  }
  ""
}

# ---------- 工具: 路径归一化 (WSL 感知) ----------
# Windows 'G:\x\y' -> (unix) '/mnt/g/x/y' / (windows) 'G:/x/y'
.env_norm <- function(p) {
  if (!nzchar(p)) return("")
  p <- gsub("\\\\", "/", p)
  if (.Platform$OS.type == "unix" && grepl("^[A-Za-z]:/", p)) {
    p <- paste0("/mnt/", tolower(substr(p, 1, 1)), "/", sub("^[A-Za-z]:/", "", p))
  }
  p
}

# ---------- 定位项目根 (含 .env 的目录) ----------
.env_find_root <- function() {
  # 优先从 --file 定位主脚本目录, 逐级向上找 .env
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArg <- sub("^--file=", "", cmdArgs[grepl("^--file=", cmdArgs)])
  start <- if (length(fileArg) && nzchar(fileArg[1])) {
    normalizePath(dirname(fileArg[1]), winslash = "/")
  } else {
    normalizePath(getwd(), winslash = "/")
  }
  cur <- start
  for (i in 1:6) {
    if (file.exists(file.path(cur, ".env"))) return(cur)
    parent <- dirname(cur)
    if (parent == cur) break
    cur <- parent
  }
  # 兜底: 从 getwd 向上
  cur <- normalizePath(getwd(), winslash = "/")
  for (i in 1:6) {
    if (file.exists(file.path(cur, ".env"))) return(cur)
    parent <- dirname(cur)
    if (parent == cur) break
    cur <- parent
  }
  stop("_env.R: cannot locate project .env (searched from script dir and getwd)")
}

.PT_PROJECT_ROOT <- .env_find_root()

# ---------- 导出核心路径 ----------
DATASET_PATH      <- .env_norm(.env_parse("DATASET_PATH"))
DATASET_PATH_TEMP <- .env_norm(.env_parse("DATASET_PATH_TEMP"))

# 存在性断言：宁可 stop() 也不要在错误路径下 dir.create(recursive=TRUE),
# 否则会在数据盘上凭空建出一棵假目录树, 把结果写进去而真实数据从未被读到。
if (!nzchar(DATASET_PATH)) {
  stop(sprintf(
    paste0("\n============================================================\n",
           " _env.R: DATASET_PATH is not set in %s\n",
           "------------------------------------------------------------\n",
           " Fix: set DATASET_PATH in that .env file.\n",
           " Refusing to fall back silently.\n",
           "============================================================\n"),
    file.path(.PT_PROJECT_ROOT, ".env")))
}
if (!dir.exists(DATASET_PATH)) {
  stop(sprintf(
    paste0("\n============================================================\n",
           " _env.R: DATASET_PATH does not exist\n",
           "------------------------------------------------------------\n",
           " .env file      : %s\n",
           " resolved value : %s\n",
           "------------------------------------------------------------\n",
           " Fix: edit DATASET_PATH in the .env file above so that it points\n",
           "      at an existing dataset root.\n",
           " Refusing to continue: the analysis scripts call\n",
           " dir.create(recursive = TRUE) on these paths, which would create a\n",
           " bogus directory tree and write results where nobody will see them.\n",
           "============================================================\n"),
    file.path(.PT_PROJECT_ROOT, ".env"), DATASET_PATH))
}

RAW       <- file.path(DATASET_PATH, "raw")
PROCESSED <- file.path(DATASET_PATH, "data/processed")
TABLES    <- file.path(DATASET_PATH, "results/tables")
FIGURES   <- file.path(DATASET_PATH, "results/figures")
TEMP      <- if (nzchar(DATASET_PATH_TEMP)) DATASET_PATH_TEMP else DATASET_PATH

# 临时区可按需创建 (不作为致命错误)
ensure_temp <- function() {
  if (nzchar(TEMP) && !dir.exists(TEMP)) dir.create(TEMP, recursive = TRUE)
  TEMP
}

invisible(NULL)
