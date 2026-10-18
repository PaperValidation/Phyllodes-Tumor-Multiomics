# 分析侧数字口径表（写作交接用）

**项目：** 乳腺叶状肿瘤多组学
**生成日期：** 2026-10-18
**数据快照：** 2026-07-30（`results/tables` 既有 15 表 mtime）
**用途：** 本表是正文/图表**唯一允许使用的数字来源**。凡本表未列的数字，不得写入稿件。
**状态：** 两处阻断项已由用户拍板（趋势=19 主/1 敏感性；可药靶=规则 A），分析侧已全部重算完毕。

---

## 一、队列与样本量（可直接写）

| 数据集 | 平台 | n | 分组 | 来源 |
|---|---|---|---|---|
| GSE231574 | EPIC 850K | 33 | Malignant 18 / Borderline 9 / Benign 2 / FA 2 / Metaplastic 2 | pheno.rds |
| GSE179458 | 450K | 96 | PT 子集 33（Malignant 10 / Benign 15 / Borderline 8） | pheno.rds |
| GSE235267 | Agilent lncRNA+mRNA | 16 | FA 4 / Benign 6 / Malignant 6 | series matrix |
| GSE78071 | Illumina HT-12 | 48 | FA 24（Simple 19 + Complex 5）/ Benign 16 / Borderline 4 / Malignant 4 | series matrix |

**探针数：** EPIC 519,767；450K 361,964。

---

## 二、差异表达（GSE78071，FDR<0.05，limma + BH）

| 指标 | 数值 | 来源文件 | 备注 |
|---|---|---|---|
| PT vs FA | **737 探针**（= 674 唯一基因） | `gse78071_DE_PT_vs_FA.csv` | 上/下调：281 / 456 探针 |
| Mal vs Ben | **62 探针**（= 62 基因） | `gse78071_DE_Mal_vs_Ben.csv` | 上调 7 / 下调 55 |
| **Ordinal trend（连续 severity）** | **19 探针（19 基因）** | `gse78071_DE_ordinal_trend_pooledvar.csv` | 主口径；17 下调 / 2 上调 |
| Ordinal trend（PT-only 敏感性） | **1 探针（ABO）** | `gse78071_DE_ordinal_trend.csv` | n=24，df 25.5，仅敏感性分析 |

**必须遵守的表述纪律：**
1. **19 是 62 的子集**，正文必须写明 "all nested within the 62 probes of the malignant-vs-benign contrast"，否则读者会误算成 81 个独立基因。
2. 旧的 `gse78071_DE_Trend.csv` **不可使用**：其对比为 `(Bor−Ben)+(Mal−Bor) ≡ Mal−Ben`，与 Mal_vs_Ben 逐字节相同（MD5 `bebef730…`）。
3. 737 是**探针**数，写"基因"须改为 674。
4. 不得把 Mal vs Ben（Mal n=4，名义富集 1.21×）写成"恶性特异表达程序"的确定结论，只能作 exploratory。

---

## 三、差异甲基化（FDR<0.05，limma on M-values）

| 对比 | 显著 DMP | 来源 |
|---|---|---|
| EPIC Bor vs Ben | **0** | `gse231574_DMP_Bor_vs_Ben.csv` |
| EPIC Mal vs Ben | **0** | `gse231574_DMP_Mal_vs_Ben.csv` |
| EPIC Mal vs Bor | **1**（min adj.P = 0.0261） | `gse231574_DMP_Mal_vs_Bor.csv` |
| EPIC ordinal trend | **0**（min adj.P = 0.7426） | `gse231574_DMP_ordinal_trend.csv` |
| EPIC DMR（DMRcate） | **0** | `gse231574_DMR_Mal_vs_Ben.csv`（仅表头） |
| 450K Mal vs Ben | **0**（min adj.P = 0.1399） | `gse179458_DMP_Mal_vs_Ben.csv` |

**表述纪律：**
- 这是**零结果**，手稿原本已如实报告，正确。必须保留配套的效能分析（EPIC benign n=2 是主因）。
- 旧 `gse231574_DMP_Trend.csv` **不可使用**（同为 Mal_vs_Ben 副本，MD5 `14596fcd…`）。
- **不要**用"名义 p<0.05 富集低于 1"来论证无信号（低 df 造成的校准假象）；用 FDR 计数。

---

## 四、甲基化–表达一致性（concordance）

| 指标 | 数值 | 来源 |
|---|---|---|
| overlap 基因数 | **13,389** | `audit_stats.py` 输出 |
| observed coherent | **51.2%**（0.5124） | 同上 |
| permutation null mean | **50.1%**（0.5013），95% 区间 0.4929–0.5097 | 同上 |
| permutation p | **0.0045**（种子 20260730；5 种子区间 0.0015–0.0080） | 同上 |
| sign correlation | −0.0225 | 同上 |

**口径说明：** 13,389 = `gse78071_meth_module_concordance.csv` 的 18,544 个 symbol ∩ 两甲基化队列 DMP 基因。**不是** 18,544。表述须保留 "weakly but significantly above chance (Δ=1.1pp)"，不可写成强共调控。

---

## 五、lncRNA-mRNA 网络（GSE235267）

| 指标 | 数值 | 来源 |
|---|---|---|
| 进展相关 lncRNA | 21（FDR<0.05 仅 3 个 + \|logFC\| top-20 放宽） | `gse235267_lncRNA_mRNA_network.csv` |
| 边数 | 1,050（其中 666 条 \|ρ\|>0.8） | 同上 |
| GO 术语 | 1,253 | `gse235267_lncRNA_partner_GO.csv` |
| GSE235267 DE（FDR<0.05） | PT vs FA 374；Mal vs Ben 103 | `gse235267_DE_*.csv` |

**纪律：** 21 个 lncRNA 的种子选择方法必须如实披露（FDR 仅 3 个，其余靠 |logFC| top-20 放宽）；网络无全基因组 FDR 校正，须在 Limitations 声明。

---

## 六、可药靶（✅ 已改为可复现规则 A）

**规则：** `adjP_PT_vs_FA < 0.05 AND pt_rel AND n_drugs > 0`
**来源：** `gse_pt_druggable_selected8_ruleA.csv`（脚本 `scripts/python/09_select_druggable_ruleA.py`）
**候选池：** 851 行（817 真基因 + 34 个 Agilent 探针 ID）

| # | 基因 | 机制轴 | log2FC | adj.P (PTvsFA) | 候选药物 | 药物状态 |
|---|---|---|---|---|---|---|
| 1 | TGFB2 | TGFb_ECM | +1.592 | 0.0174 | galunisertib; vactosertib | investigational |
| 2 | PDGFD | Stromal_other | −0.924 | 0.0192 | "PDGFR inhibitors(context)" | **无具名药物** |
| 3 | TWIST1 | Stromal_other | +0.795 | 0.0251 | "investigational(EMT)" | **无具名药物** |
| 4 | CDKN2A | CellCycle_CDK | +1.279 | 0.0385 | "CDK4/6 inhibitors(context)" | **无具名药物** |
| 5 | CCNE1 | CellCycle_CDK | +1.051 | 0.0387 | "CDK4/6 inhibitors(investigational)" | **无具名药物** |
| 6 | KIT | RTK_PDGFR_FGFR | −0.876 | 0.0421 | imatinib; sunitinib; avapritinib | **3 个均获批** |
| 7 | TGFBR2 | TGFb_ECM | −0.719 | 0.0429 | galunisertib | investigational |
| 8 | PDGFA | RTK_PDGFR_FGFR | −1.049 | 0.0493 | imatinib; nintedanib | **均获批** |

**必须遵守的表述纪律（4 条）：**
1. **DE 支持 8/8 全部 FDR<0.05**，但**全部来自 PT vs FA（实体对比）**；在 Mal vs Ben 分级对比中 **0/8 显著**。**不得**表述为随分级上调或进展信号。
2. **仅 2/8（KIT、PDGFA）有已批准的直接靶向药**；4/8 连具名化合物都没有。**"actionable agent" 一词不可用于全部 8 个**。建议表述：*"eight targets meeting a reproducible selection rule (differential expression vs fibroadenoma + PT-mechanism annotation + curated candidate agent); approved direct inhibitors exist for KIT and PDGFA only"*。
3. **删除 "independently supported"**：候选池本身取自同一 DE 表，全表 79.2% 行 DE_sig=True，属循环论证。
4. 字典条目表述：**53 gene–axis entries**（展开 111 gene–agent 对、81 唯一药名；其中 4 条无真实药物）。旧稿 "53 gene–agent–axis" 口径不准确。

**旧 8 靶点（已弃用，不得再写入正文）：** COL1A1、MMP2、PARP1 已被规则淘汰（DE 不显著，adj.P = 0.054 / 0.831 / 0.285）。若需保留 losartan/COL1A1 抗纤维化叙事，只能放 Discussion 作机制讨论，**不得列入 prioritized leads 表**。

**已修正：** 药名 `talazaparib` → **talazoparib**（源头 `07_annotate_drugs.R:58`、`rebuild_druggable.py:122` 两处均已改）。

---

## 七、配图状态（已重出，2026-10-18 19:50）

| 图 | 内容 | 变更 |
|---|---|---|
| Fig1_methylation | A: DMP 计数（519,767/361,964 探针）；B: concordance 51.2% vs 50.1%, p=0.008 | 无变化（数值本就正确） |
| Fig2_network_GO | lncRNA 网络 + GO | 无变化 |
| **Fig3_expression** | A: PT vs FA 火山图（737）；**B: 737 / 62 / 19** | ⚠️ **B 的第三根柱已从旧的 62（Trend 副本）改为 19（真 ordinal 回归）** |
| **Fig4_targets** | 规则 A 的 8 靶点；实心圆=有已批准药 | ⚠️ **8 个靶点已换（删 COL1A1/MMP2/PARP1，加 KIT/TGFBR2/PDGFA）** |

文件：`docs/manuscript/figures/*.png|pdf`，分辨率 ≈2500–3600 × 1300–1470 px（满足双栏 300 dpi）。

---

## 八、仍待用户提供的占位符

| 位置 | 内容 |
|---|---|
| §8 Data Availability | 仓库 URL（项目目前不是 git 仓库） |
| §9 CRediT | 13 个角色真名 |
| 标题页 | 作者 / 单位 / 通讯邮箱 / ORCID |

---

## 九、复现命令

```bash
# 趋势重跑（真 ordinal 回归）
MSYS_NO_PATHCONV=1 wsl.exe -d Ubuntu-24.04 -e bash -lc \
  'Rscript "/mnt/c/Users/tidus/PycharmProjects/Dermatology-Phyllodes-tumor/scripts/R/_ordinal_trend.R"'

# 可药靶规则 A 筛选
python scripts/python/09_select_druggable_ruleA.py

# concordance 置换检验
python scripts/python/audit_stats.py

# 重出全部配图
python scripts/python/08_figures.py
```

---

*本表由分析侧产出，不含正文文本。正文撰写请交由写作类专家团，并以本表数字为准。*
