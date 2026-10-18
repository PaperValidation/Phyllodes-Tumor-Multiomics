import csv, os, random, numpy as np
random.seed(20260730); np.random.seed(20260730)
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from _env import TABLES as TAB, DATASET_PATH

# ---------- 1) signature denominator ----------
# NOTE: the malignant-signature/leakage analysis was correctly abandoned at
# n=4 malignant samples (see manuscript Limitations), so the CSV no longer
# exists in TABLES. Keep this block so the script still runs end-to-end.
p = os.path.join(TAB, "gse78071_malignant_signature_score.csv")
sig = []
if os.path.exists(p):
    with open(p, newline='') as f:
        r = csv.DictReader(f)
        sig = list(r)
    print("=== SIGNATURE ===")
    print("columns:", list(sig[0].keys()))
    from collections import Counter
    cls = Counter(row.get("class", "?") for row in sig)
    print("class counts:", dict(cls))
    # malignant median
    mal = [float(row["score"]) for row in sig if row.get("class") == "Malignant_PT"]
    malmed = float(np.median(mal)) if mal else float("nan")
    print("n Malignant_PT:", len(mal), "median score:", round(malmed, 3))
    hot = [row for row in sig if row.get("class") in ("Benign_PT", "Borderline_PT") and float(row["score"]) > malmed]
    print("benign/borderline samples:", sum(1 for row in sig if row.get("class") in ("Benign_PT", "Borderline_PT")),
          "| malignant-like leakage:", len(hot))
    if hot:
        print("  leaked samples:", [(h["sample"], h["class"]) for h in hot])
else:
    print("=== SIGNATURE === (skipped: gse78071_malignant_signature_score.csv not found; "
          "leakage analysis abandoned at n=4 malignant)")

# ---------- 2) concordance permutation null ----------
# methylation direction: UNION of GSE231574 + GSE179458 DMP Mal_vs_Ben (logFC = methylation, Mal - Ben)
# expression direction: GSE78071 PT_vs_FA (logFC = expression, PT - FA) from concordance CSV
def load_symbol_logfc(path, genecol, lfccol):
    d = {}
    with open(path, newline='') as f:
        for row in csv.DictReader(f):
            g = row.get(genecol, "")
            if g in ("", "NA", None):
                continue
            try:
                v = float(row[lfccol])
            except:
                continue
            d.setdefault(g, []).append(v)
    return {g: float(np.mean(vs)) for g, vs in d.items()}

m1 = load_symbol_logfc(os.path.join(TAB, "gse231574_DMP_Mal_vs_Ben.csv"), "gene", "logFC")
m2 = load_symbol_logfc(os.path.join(TAB, "gse179458_DMP_Mal_vs_Ben.csv"), "gene", "logFC")
meth = {}
for g in set(m1) | set(m2):
    vs = []
    if g in m1: vs.append(m1[g])
    if g in m2: vs.append(m2[g])
    meth[g] = float(np.mean(vs))
expr = load_symbol_logfc(os.path.join(TAB, "gse78071_meth_module_concordance.csv"), "symbol", "logFC")
genes = sorted(set(meth) & set(expr))
print("\n=== CONCORDANCE PERMUTATION (union of both methylation cohorts) ===")
print("overlap genes (meth-direction x expr-direction):", len(genes))
mvals = np.array([meth[g] for g in genes])
evalsv = np.array([expr[g] for g in genes])

def consistency(m, e, coherent=True):
    # coherent: hypermethylation (m>0) suppresses expression (e<0)
    if coherent:
        agree = ((m > 0) & (e < 0)) | ((m < 0) & (e > 0))
    else:
        agree = ((m > 0) & (e > 0)) | ((m < 0) & (e < 0))
    return agree.mean()

obs_coh = consistency(mvals, evalsv, True)
obs_same = consistency(mvals, evalsv, False)
print("observed coherent-consistency (hyper<->down):", round(obs_coh, 4))
print("observed same-sign consistency:", round(obs_same, 4))

# permutation: shuffle expression signs across genes (keep methylation fixed)
N = 2000
null_coh = np.empty(N)
es = evalsv.copy()
for i in range(N):
    perm = np.random.permutation(evalsv)
    null_coh[i] = consistency(mvals, perm, True)
null_mean = null_coh.mean()
pval = (np.sum(null_coh >= obs_coh) + 1) / (N + 1)
print(f"permutation null mean (coherent): {null_mean:.4f}  [95% range {np.percentile(null_coh,2.5):.4f}-{np.percentile(null_coh,97.5):.4f}]")
print(f"observed vs null: p(perm>=obs) = {pval:.4f}")
# also correlation between meth and expr direction (existing global anticorrelation)
rho = np.corrcoef(np.sign(mvals), np.sign(evalsv))[0, 1]
print("sign-correlation(meth,expr) =", round(rho, 4))

# reconstruct the pipeline 72.9% style: how many genes had BOTH a directional meth call
# (already all do here). Compare to reported 18544/25434.
print("\noverlap genes used for concordance:", len(genes))
