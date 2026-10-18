"""Single source of truth for project data paths.

Reads DATASET_PATH / DATASET_PATH_TEMP from the repo `.env` (the same
variables the R pipeline and run_pipeline.sh already honour), so every
script points at the same dataset root instead of a hard-coded local copy.

Resolves DATASET_PATH into a path usable from both Windows Python (`G:/...`)
and WSL Python (`/mnt/g/...`).

Exposes:
    DATASET_PATH      - root dataset dir ( == .env DATASET_PATH )
    DATASET_PATH_TEMP - temp dataset dir ( == .env DATASET_PATH_TEMP )
    RAW               - DATASET_PATH/raw
    PROCESSED         - DATASET_PATH/data/processed
    TABLES            - DATASET_PATH/results/tables
    FIGURES           - DATASET_PATH/results/figures
    TEMP              - DATASET_PATH_TEMP (temp/intermediate scratch)
"""
import os
import re
import sys

_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _read_dotenv():
    env = {}
    p = os.path.join(_PROJECT_ROOT, ".env")
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                k, v = line.split("=", 1)
                env[k.strip()] = v.strip().strip('"').strip("'")
    return env


def _norm(p):
    """Normalise a Windows/WSL path for the current platform."""
    if not p:
        return ""
    p = p.replace("\\", "/").strip()
    if sys.platform.startswith("linux"):
        m = re.match(r"^([A-Za-z]):/(.*)$", p)
        if m:
            p = "/mnt/" + m.group(1).lower() + "/" + m.group(2)
    return os.path.normpath(p) if p else ""


def _env_file():
    return os.path.join(_PROJECT_ROOT, ".env")


def _die(key, raw, resolved, reason):
    """Fail loudly instead of silently falling back to the code repo.

    Silent fallback is the worst failure mode: results get written into the
    project repository (or a freshly invented directory tree) and nobody
    notices until the numbers are wrong.
    """
    raise RuntimeError(
        "\n"
        "============================================================\n"
        " _env.py: cannot resolve %s\n"
        "------------------------------------------------------------\n"
        " .env file      : %s\n"
        " resolved value : %r\n"
        " problem        : %s\n"
        "------------------------------------------------------------\n"
        " Fix: edit DATASET_PATH / DATASET_PATH_TEMP in the .env file above\n"
        "      so that it points at an existing dataset root, or export the\n"
        "      variable in the shell before running the script.\n"
        " Refusing to fall back silently to the source tree.\n"
        "============================================================\n"
        % (key, _env_file(), raw if raw else resolved, reason)
    )


def dataset_path():
    # 1) shell-exported env var (WSL / run_pipeline.sh export)
    raw = os.environ.get("DATASET_PATH", "").strip()
    src = "environment variable DATASET_PATH"
    # 2) parsed from .env
    if not raw:
        raw = _read_dotenv().get("DATASET_PATH", "").strip()
        src = ".env DATASET_PATH"
    if not raw:
        _die("DATASET_PATH", raw, raw,
             "not set in %s and not exported in the environment" % _env_file())
    p = _norm(raw)
    if not os.path.isdir(p):
        _die("DATASET_PATH", raw, p,
             "directory does not exist (from %s)" % src)
    return p


def dataset_path_temp():
    # 1) shell-exported env var
    raw = os.environ.get("DATASET_PATH_TEMP", "").strip()
    # 2) parsed from .env
    if not raw:
        raw = _read_dotenv().get("DATASET_PATH_TEMP", "").strip()
    # 3) not configured -> reuse the dataset root (documented behaviour)
    if not raw:
        return DATASET_PATH
    p = _norm(raw)
    # Temp area may legitimately not exist yet; it is created on demand by
    # ensure_temp(). Do NOT silently swap it for some other directory.
    return p


def ensure_temp():
    """Create DATASET_PATH_TEMP if missing and return it."""
    p = dataset_path_temp()
    if p and not os.path.isdir(p):
        os.makedirs(p, exist_ok=True)
    return p


DATASET_PATH = dataset_path()
DATASET_PATH_TEMP = dataset_path_temp()
RAW = os.path.join(DATASET_PATH, "raw")
PROCESSED = os.path.join(DATASET_PATH, "data", "processed")
TABLES = os.path.join(DATASET_PATH, "results", "tables")
FIGURES = os.path.join(DATASET_PATH, "results", "figures")
TEMP = DATASET_PATH_TEMP
