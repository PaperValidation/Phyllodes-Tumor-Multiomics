#!/usr/bin/env python
# ============================================================
# 09_select_druggable_ruleA.py
#   可药靶优先化的【可复现规则】版本（方案 A，2026-10-18 用户拍板采纳）
#
# 背景
#   旧稿的 8 个优先靶点来自 rebuild_druggable.py 中手写的 MECH_VALID 集合
#   （{MMP2,TGFB2,TWIST1,CCNE1,CDKN2A,COL1A1,PDGFD} + PARP1），属人工圈定，
#   但正文写成 "eight genes received at least one actionable agent" —— 按该
#   描述实测应得 53 个，属表述与数据不符。审稿人三行代码即可发现。
#
# 本脚本用**完全可复现的显式规则**重新筛选：
#       DE_sig_PT_vs_FA (adj.P.Val < 0.05)
#   AND pt_rel            (有 PT 机制轴赋值)
#   AND n_drugs > 0       (curated 字典中有 >=1 个候选化合物)
#
# 输出（新文件，绝不覆盖既有表）：
#   gse_pt_druggable_selected8_ruleA.csv   规则选出的靶点 + 完整注释
#   _druggable_ruleA_summary.txt           运行摘要
#
# 只读既有 candidate 表；不重跑任何上游分析。
# ============================================================
import csv
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from _env import TABLES as TAB  # noqa: E402

SRC = os.path.join(TAB, "gse_pt_druggable_candidates.csv")
OUT = os.path.join(TAB, "gse_pt_druggable_selected8_ruleA.csv")
LOG = os.path.join(TAB, "_druggable_ruleA_summary.txt")

# ------------------------------------------------------------
# 药物获批状态（curated；用于诚实标注 "actionable" 一词是否成立）
# 只收录本规则会触及的药物名，避免过度推断。
# approved       = 已有监管批准（可能为其他适应症）
# investigational= 未获批 / 试验中止
# ------------------------------------------------------------
AGENT_STATUS = {
    # 已获批
    "imatinib": "approved", "sunitinib": "approved", "avapritinib": "approved",
    "nintedanib": "approved", "pazopanib": "approved", "regorafenib": "approved",
    "axitinib": "approved", "olaparib": "approved", "niraparib": "approved",
    "talazoparib": "approved", "losartan": "approved(non-oncology)",
    "pirfenidone": "approved(antifibrotic)", "doxycycline": "approved(non-oncology)",
    "everolimus": "approved", "palbociclib": "approved", "ribociclib": "approved",
    "abemaciclib": "approved",
    # 未获批
    "galunisertib": "investigational", "vactosertib": "investigational",
    "marimastat": "investigational", "prinomastat": "investigational",
    "fresolimumab": "investigational", "pamrevlumab": "investigational",
    "simtuzumab": "investigational", "dinaciclib": "investigational",
    "seliciclib": "investigational", "capivasertib": "investigational",
    "ipatasertib": "investigational",
}


def agent_status(drugs_field):
    """Return (n_approved, n_invest, status_label) for a 'a;b;c' drug string.

    Strings like 'CDK4/6 inhibitors(context)' / 'investigational(EMT)' / 'no approved'
    contain no named agent -> classified as 'no-named-agent'.
    """
    names = []
    for part in (drugs_field or "").split(";"):
        p = part.strip()
        if not p:
            continue
        # strip parenthetical qualifiers: 'marimastat(investigational)' -> 'marimastat'
        base = p.split("(")[0].strip().lower()
        names.append((p, base))
    if not names:
        return 0, 0, "none"
    approved = sum(1 for _, b in names if AGENT_STATUS.get(b, "").startswith("approved"))
    invest = sum(1 for _, b in names if AGENT_STATUS.get(b, "") == "investigational")
    named = sum(1 for _, b in names if b in AGENT_STATUS)
    if named == 0:
        return 0, 0, "no-named-agent"
    if approved > 0:
        return approved, invest, "has-approved"
    if invest > 0:
        return 0, invest, "investigational-only"
    return 0, 0, "unclear"


def main():
    if not os.path.exists(SRC):
        raise SystemExit(f"missing candidate table: {SRC}")

    with open(SRC, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    print(f"read {len(rows)} candidate rows from {os.path.basename(SRC)}")

    # ---- apply the reproducible rule ----
    def fnum(v, default=1.0):
        try:
            return float(v)
        except (TypeError, ValueError):
            return default

    sel = []
    for r in rows:
        adjp = fnum(r.get("adjP_PT_vs_FA"), 1.0)
        pt_rel = (r.get("pt_rel", "") or "").strip().lower() == "true"
        try:
            ndrug = int(r.get("n_drugs") or 0)
        except ValueError:
            ndrug = 0
        if adjp < 0.05 and pt_rel and ndrug > 0:
            sel.append(r)

    # ---- annotate ----
    out_rows = []
    for r in sorted(sel, key=lambda x: fnum(x.get("adjP_PT_vs_FA"))):
        n_app, n_inv, label = agent_status(r.get("drugs", ""))
        rr = dict(r)
        rr["n_agents_approved"] = n_app
        rr["n_agents_investigational"] = n_inv
        rr["agent_status"] = label
        out_rows.append(rr)

    fields = list(out_rows[0].keys()) if out_rows else []
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)

    # ---- summary ----
    lines = []
    lines.append("=" * 69)
    lines.append(" Druggable target selection — REPRODUCIBLE RULE (option A)")
    lines.append("   rule: adjP_PT_vs_FA < 0.05  AND  pt_rel  AND  n_drugs > 0")
    lines.append("=" * 69)
    lines.append(f"source candidates : {len(rows)}")
    n_drug_any = sum(1 for r in rows if int(r.get('n_drugs') or 0) > 0)
    n_ptrel = sum(1 for r in rows if (r.get('pt_rel', '') or '').strip().lower() == 'true')
    n_sig = sum(1 for r in rows if fnum(r.get('adjP_PT_vs_FA')) < 0.05)
    lines.append(f"  with >=1 curated agent : {n_drug_any}")
    lines.append(f"  with PT-mechanism axis : {n_ptrel}")
    lines.append(f"  DE-significant PTvsFA  : {n_sig}")
    lines.append(f"SELECTED (all 3)         : {len(out_rows)}")
    lines.append("")
    hdr = f"{'gene':<10}{'axis':<18}{'log2FC':>9}{'adjP':>11}{'n_drug':>7}  {'status':<20}drugs"
    lines.append(hdr)
    lines.append("-" * 108)
    for r in out_rows:
        lines.append(
            f"{r['gene']:<10}{r.get('pt_pathway',''):<18}"
            f"{fnum(r.get('logFC_PT_vs_FA'), 0):>9.3f}{fnum(r.get('adjP_PT_vs_FA')):>11.4f}"
            f"{r['n_drugs']:>7}  {r['agent_status']:<20}{r.get('drugs','')[:40]}"
        )
    lines.append("")
    n_app = sum(1 for r in out_rows if r["agent_status"] == "has-approved")
    lines.append(f"targets with >=1 approved agent : {n_app}/{len(out_rows)}")
    lines.append(f"targets with only investigational/context agents : {len(out_rows)-n_app}/{len(out_rows)}")
    lines.append("")
    lines.append("CAVEAT for the manuscript:")
    lines.append("  All DE support derives from the PT-vs-FA (entity-level) contrast.")
    for r in out_rows:
        pass
    n_mb_sig = sum(1 for r in out_rows if fnum(r.get("adjP_Mal_vs_Ben"), 1.0) < 0.05)
    lines.append(f"  Targets also significant in Mal-vs-Ben (grade contrast): {n_mb_sig}/{len(out_rows)}")
    lines.append("")
    lines.append("WROTE: " + OUT)
    txt = "\n".join(lines) + "\n"
    with open(LOG, "w", encoding="utf-8") as f:
        f.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
