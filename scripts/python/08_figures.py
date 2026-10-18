#!/usr/bin/env python
# ============================================================
# 08_figures.py — 生成 manuscript 的 Fig1-4 + Table 支撑图 (Python 版)
# 数据来源: DATASET_PATH/results/tables/*.csv (由 scripts/_env.py 从 .env 解析, 与 R 流水线 02-07 同源)
# 依赖: matplotlib, networkx, pandas, numpy, scipy
# 说明: 原稿计划用 R/08_figures.R; 本环境无 R, 用 Python 生成等价图文件,
#       数值与稿件 manuscript_v1.md 一致 (已 verify_pt.py 校验).
# ============================================================
import os, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import sys, os
_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.dirname(_HERE))
from _env import TABLES as TAB, DATASET_PATH
# 手稿配图输出目录: 随项目走, 与机器/用户名无关
FIG = os.path.join(_PROJECT_ROOT, "docs", "manuscript", "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 10, "figure.dpi": 300, "axes.grid": False,
                     "font.family": "DejaVu Sans"})

def load(p):
    with open(p, newline='') as f:
        return list(csv.DictReader(f))

# ---------------- Fig 1: Methylation (DMP near-null + concordance) ----------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
contrasts = ["EPIC\nBor_vs_Ben", "EPIC\nMal_vs_Ben", "EPIC\nMal_vs_Bor",
             "EPIC\nTrend", "450k\nMal_vs_Ben"]
dmp_counts = [0, 0, 1, 0, 0]
b = ax[0].bar(contrasts, dmp_counts, color="#4C72B0", edgecolor="black", width=0.6)
ax[0].set_ylabel("Significant DMPs (FDR<0.05)")
ax[0].set_title("A. CpG-level differential methylation is near-null\n(519,767 EPIC / 361,964 450k probes tested)")
ax[0].set_ylim(0, 3)
for i, v in enumerate(dmp_counts):
    ax[0].text(i, v + 0.05, str(v), ha="center", fontsize=9)
ax[0].text(0.5, -0.28, "EPIC benign n=2 → Mal_vs_Ben / Bor_vs_Ben contrasts underpowered (see power analysis)",
           transform=ax[0].transAxes, ha="center", fontsize=7.5, color="#8B0000")

# concordance permutation bar
obs, nullv = 51.2, 50.1
ax[1].bar(["Observed\n(13,389 genes)", "Permutation\nnull (2,000x)"], [obs, nullv],
          color=["#4C72B0", "#BBBBBB"], edgecolor="black", width=0.55)
ax[1].set_ylim(0, 70)
ax[1].set_ylabel("Directional concordance (%)")
ax[1].set_title("B. Methylation-expression directional concordance\nweakly above chance (51.2% vs 50.1%, p=0.008)")
for i, v in enumerate([obs, nullv]):
    ax[1].text(i, v + 1, f"{v}%", ha="center", fontsize=9)
ax[1].text(0.5, -0.30, "Cross-platform direction test; small but significant effect (p<0.01, \u0394=1.1pp)",
           transform=ax[1].transAxes, ha="center", fontsize=7.5)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "Fig1_methylation.png"), bbox_inches="tight")
plt.savefig(os.path.join(FIG, "Fig1_methylation.pdf"), bbox_inches="tight")
plt.close()

# ---------------- Fig 2: lncRNA network + GO ----------------
import networkx as nx
net = load(os.path.join(TAB, "gse235267_lncRNA_mRNA_network.csv"))
edges = [(r["lncRNA"], r["mRNA"], abs(float(r["rho"]))) for r in net]
sub = [(u, v, w) for (u, v, w) in edges if w > 0.8]   # 666 strongest edges
G = nx.Graph()
for u, v, w in sub:
    G.add_edge(u, v, weight=w)
lnc_nodes = set(r["lncRNA"] for r in net)
pos = nx.spring_layout(G, k=0.25, iterations=60, seed=7)
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
node_colors = ["#C44E52" if n in lnc_nodes else "#4C72B0" for n in G.nodes()]
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=28, ax=ax[0], alpha=0.85)
ew = [G[u][v]["weight"] for u, v in G.edges()]
nx.draw_networkx_edges(G, pos, edge_color=ew, edge_cmap=plt.cm.viridis,
                       width=0.4, alpha=0.5, ax=ax[0])
ax[0].set_title(f"A. lncRNA-mRNA co-expression network\n(strongest 666 edges, |rho|>0.8 of 1,050 total; 21 lncRNAs)")
ax[0].axis("off")
ax[0].legend(handles=[Patch(color="#C44E52", label="lncRNA (21)"),
                      Patch(color="#4C72B0", label="mRNA partner")], loc="upper right", fontsize=8)

go = load(os.path.join(TAB, "gse235267_lncRNA_partner_GO.csv"))
go_sorted = sorted(go, key=lambda r: float(r["p.adjust"]))[:10]
terms = [f"{r['Term']}" for r in go_sorted][::-1]
padj = [-math.log10(float(r["p.adjust"])) for r in go_sorted][::-1]
ax[1].barh(range(len(terms)), padj, color="#8172B3", edgecolor="black")
ax[1].set_yticks(range(len(terms)))
ax[1].set_yticklabels([t[:34] for t in terms], fontsize=7.5)
ax[1].set_xlabel("-log10(adjusted p)")
ax[1].set_title("B. Top GO terms (partner mRNAs)\n(1,253 enriched terms total)")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "Fig2_network_GO.png"), bbox_inches="tight")
plt.savefig(os.path.join(FIG, "Fig2_network_GO.pdf"), bbox_inches="tight")
plt.close()

# ---------------- Fig 3: Expression validation (volcano + DE counts) ----------------
de = load(os.path.join(TAB, "gse78071_DE_PT_vs_FA.csv"))
x = np.array([float(r["logFC"]) for r in de if r["logFC"] not in ("", "NA")])
y = np.array([min(-math.log10(float(r["adj.P.Val"])), 300) for r in de if r["adj.P.Val"] not in ("", "NA")])
sig = np.array([float(r["adj.P.Val"]) < 0.05 for r in de if r["adj.P.Val"] not in ("", "NA")])
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
ax[0].scatter(x[~sig], y[~sig], s=4, color="#BBBBBB", alpha=0.4, label="n.s.")
ax[0].scatter(x[sig], y[sig], s=6, color="#C44E52", alpha=0.7, label="FDR<0.05")
ax[0].axhline(-math.log10(0.05), color="black", ls="--", lw=0.8)
ax[0].set_xlabel("log2(PT/FA) expression fold-change\n(limma, log2 + quantile normalized)")
ax[0].set_ylabel("-log10(adjusted p)")
ax[0].set_title(f"A. GSE78071 PT vs FA volcano\n{int(sig.sum())} significant genes (FDR<0.05)")
ax[0].legend(fontsize=8, loc="upper right")
xm = max(abs(x.min()), abs(x.max())) * 1.05
ax[0].set_xlim(-xm, xm)

def sigcount(f):
    p = os.path.join(TAB, f)
    if not os.path.exists(p):
        return None
    return sum(1 for r in load(p)
               if r.get("adj.P.Val") not in ("", "NA") and float(r["adj.P.Val"]) < 0.05)

# NOTE (2026-10-18): the legacy "gse78071_DE_Trend.csv" is algebraically identical to
# Mal_vs_Ben -- the contrast was (Bor-Ben)+(Mal-Bor) === Mal-Ben -- so plotting it as a
# separate bar double-counted the same data. The ordinal bar now uses the genuine
# continuous-severity regression (gse78071_DE_ordinal_trend_pooledvar.csv, df 47.5,
# same variance structure as the pairwise model => directly comparable).
de_counts = [sigcount("gse78071_DE_PT_vs_FA.csv"),
             sigcount("gse78071_DE_Mal_vs_Ben.csv"),
             sigcount("gse78071_DE_ordinal_trend_pooledvar.csv")]
labels = ["PT vs FA", "Mal vs Ben", "Ordinal trend\n(continuous severity)"]
if any(v is None for v in de_counts):
    raise SystemExit(f"missing DE table for Fig3B: {de_counts}")
bb = ax[1].bar(labels, de_counts, color=["#C44E52", "#4C72B0", "#55A868"], edgecolor="black")
ax[1].set_ylabel("Significant DE genes (FDR<0.05)")
ax[1].set_title("B. Differential expression in GSE78071\n(48 samples: FA + benign/borderline/malignant PT)")
ax[1].set_ylim(0, max(de_counts) * 1.18)
for i, v in enumerate(de_counts):
    ax[1].text(i, v + max(de_counts) * 0.03, str(v), ha="center", fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "Fig3_expression.png"), bbox_inches="tight")
plt.savefig(os.path.join(FIG, "Fig3_expression.pdf"), bbox_inches="tight")
plt.close()

# ---------------- Fig 4: Druggable targets dot plot ----------------
# NOTE (2026-10-18): targets are now derived from the REPRODUCIBLE selection rule
# (scripts/python/09_select_druggable_ruleA.py):
#   adjP_PT_vs_FA < 0.05  AND  pt_rel  AND  n_drugs > 0
# This replaces the previously hand-curated MECH_VALID list, which the manuscript
# described as an "actionable agent" filter although it was manually nominated.
SEL = os.path.join(TAB, "gse_pt_druggable_selected8_ruleA.csv")
if not os.path.exists(SEL):
    raise SystemExit("missing %s -- run scripts/python/09_select_druggable_ruleA.py first" % SEL)
sel = load(SEL)
cand = {r["gene"]: r for r in load(os.path.join(TAB, "gse_pt_druggable_candidates.csv"))}
# order: approved-agent targets first (top), then investigational/context-only
sel_sorted = sorted(sel, key=lambda r: (r.get("agent_status") != "has-approved",
                                        float(r.get("adjP_PT_vs_FA") or 1)))
order = [r["gene"] for r in sel_sorted][::-1]
lfc = [float(cand[g]["log2FC_expr_PT_vs_FA"]) for g in order]
desig = [cand[g]["DE_sig_PT_vs_FA"] == "TRUE" for g in order]
has_appr = {r["gene"]: (r.get("agent_status") == "has-approved") for r in sel}
colors = ["#C44E52" if d else "#BBBBBB" for d in desig]
fig, ax = plt.subplots(figsize=(8.6, 5))
for i, (g, v, d) in enumerate(zip(order, lfc, desig)):
    appr = has_appr.get(g, False)
    marker = "o" if appr else "s"
    ax.scatter(v, i, s=120, color=colors[i], marker=marker,
               edgecolor="black", zorder=3,
               facecolors=colors[i] if appr else "none")
    tag = "" if appr else "  (no approved direct agent)"
    if not d:
        tag = " (n.s.)" + tag
    ax.annotate(g + tag, (v, i), xytext=(6, 0), textcoords="offset points",
                va="center", fontsize=8)
ax.axvline(0, color="black", lw=0.8)
ax.set_yticks(range(len(order))); ax.set_yticklabels([])
ax.set_xlabel("log2(PT/FA) expression fold-change")
ax.set_title("Druggable targets from a reproducible selection rule\n"
             "(all DE-significant PT vs FA; filled circle = >=1 approved direct agent)")
ax.set_ylim(-0.6, len(order) - 0.4)
ax.legend(handles=[
    Patch(color="#C44E52", label="DE-significant (PT vs FA)"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor="#C44E52",
           markeredgecolor="black", markersize=9, label="filled circle = approved agent exists"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor="white",
           markeredgecolor="black", markersize=9, label="hollow square = investigational / context only"),
], loc="lower right", fontsize=7.5)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "Fig4_targets.png"), bbox_inches="tight")
plt.savefig(os.path.join(FIG, "Fig4_targets.pdf"), bbox_inches="tight")
plt.close()

print("FIGURES WRITTEN to", FIG)
for f in ["Fig1_methylation", "Fig2_network_GO", "Fig3_expression", "Fig4_targets"]:
    print("  ", f + ".png / .pdf")
