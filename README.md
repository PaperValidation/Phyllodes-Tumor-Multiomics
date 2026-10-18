# Phyllodes tumor multi-omics re-analysis

Code, summary result tables and figure sources for the manuscript:

> **Integrative transcriptomic and DNA methylomic profiling reveals a stroma-dominated
> progression program and druggable targets in phyllodes tumors**

---

## What this repository contains

| Directory | Contents |
|---|---|
| `scripts/R/` | Analysis pipeline (steps 02–07) plus `_env.R` path resolver and `_ordinal_trend.R` |
| `scripts/python/` | Figure generation (`08_figures.py`), druggable-target selection (`09_select_druggable_ruleA.py`), concordance permutation test (`audit_stats.py`) |
| `scripts/bash/` | Data download (`01_download_data.sh`) and pipeline orchestration (`run_pipeline.sh`) |
| `results/tables/` | **Summary** result tables (differential expression, network, GO, druggability). See "Large intermediate tables" below |
| `results/summaries/` | Run logs for the ordinal-trend refit and the druggable-target selection rule |
| `docs/manuscript/figures/` | Manuscript figure sources (PNG + PDF, 300 dpi) |
| `docs/manuscript/review/` | Methodological provenance: the audit report and the numbers ledger that maps every reported value to its source file |

**What is intentionally not published here:** the manuscript text and the submission
package. This repository publishes the *analysis*, not the article; the peer-review
process is not pre-empted. The numbers ledger documents the provenance of every value
reported in the manuscript without reproducing the manuscript itself.

## Data sources

All four datasets are public and were retrieved from NCBI GEO via `GEOquery`:

| Dataset | Platform | n | Contrast(s) used |
|---|---|---|---|
| [GSE231574](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231574) | Illumina HumanMethylationEPIC | 33 | Malignant 18 / Borderline 9 / Benign 2 / FA 2 / Metaplastic 2 |
| [GSE179458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179458) | Illumina HumanMethylation450K | 96 | PT subset 33 (Malignant 10 / Benign 15 / Borderline 8) |
| [GSE235267](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235267) | Agilent lncRNA + mRNA | 16 | FA 4 / Benign 6 / Malignant 6 |
| [GSE78071](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78071) | Illumina HumanHT-12 | 48 | FA 24 / Benign 16 / Borderline 4 / Malignant 4 |

No new data were generated. All samples are de-identified public data.

## Reproducing the analysis

### 1. Environment

```bash
# R 4.3.3 / Bioconductor 3.18
#   minfi, DMRcate, limma, GEOquery, ConsensusClusterPlus,
#   org.Hs.eg.db, IlluminaHumanMethylation{EPIC,450k} annotation packages
# Python >= 3.10
#   pandas numpy scipy matplotlib networkx
```

### 2. Configure data paths

Create a `.env` file in the repository root:

```dotenv
DATASET_PATH='E:\path\to\Phyllodes tumor'
DATASET_PATH_TEMP='F:\path\to\temp'
```

All scripts resolve paths through `scripts/_env.py` (Python) or `scripts/R/_env.R` (R).
Both **abort with a clear error** if `DATASET_PATH` does not exist — they never fall
back silently to the source tree.

Layout expected under `DATASET_PATH`: `raw/`, `data/processed/`, `results/tables/`, `results/figures/`.

### 3. Download and run

```bash
bash scripts/bash/01_download_data.sh   # checkpointed GEO download
bash scripts/bash/run_pipeline.sh       # steps 02 -> 06
Rscript scripts/R/_ordinal_trend.R      # continuous ordinal severity models
python scripts/python/09_select_druggable_ruleA.py   # druggable target selection
python scripts/python/08_figures.py     # all manuscript figures
python scripts/python/audit_stats.py    # concordance permutation test
```

## Key methodological notes

**Expression normalisation (GSE78071).** The GEO series matrix contains no Illumina
negative-control probes (GPL14951 `Probe_Type` is entirely A/I/M/S, 0 control rows), so
`limma::neqc()`'s normexp background correction is not applicable. Intensities were
therefore `log2`-transformed and quantile-normalised between arrays
(`limma::normalizeBetweenArrays(method = "quantile")`).

**Ordinal severity model.** The trend analysis fits histological grade as a *continuous*
covariate (benign = 1, borderline = 2, malignant = 3) inside the same four-group linear
model used for the pairwise contrasts, so residual degrees of freedom and variance
structure are identical. An earlier contrast written as
`(borderline − benign) + (malignant − borderline)` was **discarded**: it is algebraically
identical to `malignant − benign` and would have duplicated that result.

**Druggable target selection.** Targets are selected by an explicit, fully reproducible
rule applied to all 851 candidates:

```
adj.P.Val(PT vs FA) < 0.05  AND  PT-mechanism annotation  AND  >= 1 curated candidate agent
```

This yields eight targets with no manual nomination step. Note that the
differential-expression support is **entity-level** (PT vs fibroadenoma): none of the
eight reaches significance in the malignant-vs-benign grade contrast. Approved direct
inhibitors exist for KIT and PDGFA only.

## Large intermediate tables

Probe-level DMP tables (GSE231574, 519,767 probes × 6 contrasts; GSE179458, 361,964
probes) are **not** stored in this repository: together they exceed 460 MB. They are
intermediate products that are fully regenerated by `scripts/bash/run_pipeline.sh` and
`scripts/R/_ordinal_trend.R` from the public GEO inputs listed above.

The summary tables required to verify every number reported in the manuscript are
included under `results/tables/`.

## Numbers used in the manuscript

`docs/numbers_ledger_*.md` is the authoritative record of every quantitative value
reported in the manuscript, together with the source file and the reporting constraints
for each section. If there is ever a discrepancy between the manuscript text and this
repository, the ledger and the CSV files here take precedence.

## License

Code: MIT (see `LICENSE`).
Result tables and figures: CC BY 4.0.

## Citation

See `CITATION.cff`. Please also cite the original data-generating publications for the
four GEO datasets.
