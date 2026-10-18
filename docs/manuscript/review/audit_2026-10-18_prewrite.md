# 投稿前产出物全面复核报告

**项目：** 乳腺叶状肿瘤（Phyllodes Tumor）多组学 — Integrative transcriptomic and DNA methylomic profiling
**复核日期：** 2026-10-18
**复核方式：** Biomni 生物医学专家团五线独立复核（环境运维 / 基因组学 / 药物学 / 代码加固），主理人逐项独立复验
**数据快照：** 2026-07-30（`results/tables` 既有 15 表的 mtime 一致）
**结论：** 数据与核心统计**可靠可用**；发现 **1 项分母级错误（趋势分析）**、**1 项叙事与数据不符（可药靶）**，二者均需改稿后才能开写。

---

## 一、总体判定

| 板块 | 判定 | 依据 |
|---|---|---|
| 队列描述与样本量（Table 1） | ✅ 可直接用 | 四队列全部独立复核吻合 |
| GSE78071 表达层（737 / 62 / 7 上调） | ✅ 可直接用 | 探针级精确复现 |
| 甲基化层（0/0/1/0 DMP；450k 零检出） | ✅ 可直接用 | 手稿原本已如实报告零结果 |
| 甲基化–表达一致性（51.2% vs 50.1%, p=0.0045） | ✅ 可直接用 | 重跑 `audit_stats.py` 完整复现 |
| **ordinal 趋势（原称 62 基因）** | ⛔ **必须重写** | 原对比代数上恒等于 Mal−Ben，属重复计数；正确值 19 |
| **可药靶 8 靶点** | ⚠️ **数字可用、表述必须重写** | 硬编码名单被写成筛选结果等 5 处 |
| 配图（Fig1–4） | ✅ 文件有效 | 分辨率 ≈3000×1300 px，达 300 dpi 双栏要求 |

---

## 二、阻断项 A：ordinal 趋势分析（已查清并重跑）

### 2.1 原缺陷（分母级）

`scripts/R/05_expr_gse78071.R:70` 与 `02_meth_gse231574.R:113` 的对比定义为：

```r
Trend      = (Borderline − Benign) + (Malignant − Borderline)
Mal_vs_Ben = Malignant − Benign
```

代数上 `(B−Ben)+(M−B) ≡ M−Ben`，**两者是同一个对比**。md5 校验证实：

| 文件对 | MD5 | 结果 |
|---|---|---|
| `gse78071_DE_Mal_vs_Ben.csv` / `gse78071_DE_Trend.csv` | 均 `bebef730…` | 逐字节相同 |
| `gse231574_DMP_Mal_vs_Ben.csv` / `gse231574_DMP_Trend.csv` | 均 `14596fcd…` | 逐字节相同 |

**后果：** 手稿摘要、§2.2、§3.1、§3.3、Discussion、Fig 3B 共 5 处，把"62 ordinal-trend"与"62 malignant-vs-benign"作为两个独立结果并列，实为同一份数据重复计数。

### 2.2 重跑结果（真连续回归，severity Ben=1/Bor=2/Mal=3）

新脚本 `scripts/R/_ordinal_trend.R`（不改旧脚本、不覆盖旧 CSV）：

| 模型 | 设计 | df | FDR<0.05 | p<0.001 富集 | mean p | 用途 |
|---|---|---|---|---|---|---|
| A 严格 PT-only | n=24（Ben16/Bor4/Mal4） | 25.5 | **1** | 1.43× | 0.5287（偏保守） | 敏感性分析 |
| **B 方差合并** | n=48，对比 `(0·FA−8·Ben+2·Bor+6·Mal)/14` | 47.5 | **19** | **4.63×** | 0.4948（校准好） | **正文主口径** |
| 旧 pairwise | Mal−Ben | 47.5 | 62 | — | 0.4805 | 参照 |

- 甲基化 ordinal：**0 个 DMP**（nominal 0.919×，min adj.P = 0.7426）。
- **19 ⊂ 62 且 1 ⊂ 62**（两模型显著集均为旧 62 的子集）。
- 新表 MD5 `7112b59f…` ≠ 旧表 `bebef730…`，确认是真回归而非副本。

### 2.3 ⚠️ 一条被纠正的判断（避免写错）

**"模型 A 名义富集 0.751× 低于随机"不能用来论证"无趋势"。** 校准检验表明模型 A 的 p 值分布系统性偏保守（mean p=0.5287 > 0.5；`frac(p<0.05)` 仅 0.0375 < 0.05），这是低残差自由度下 eBayes 收缩过度的产物。**同一批样本改用模型 B（df 47.5）富集即回到 1.15×，p<0.001 处达 4.63×。** 因此正文应使用 **FDR 计数 + 尾部富集**，不使用 0.05 名义富集倍数。

### 2.4 19 个基因的生物学方向（17 下调 / 2 上调）

| 方向 | 基因 |
|---|---|
| 下调（17） | ABO, LOC389816, MAPT, LRRC26, **PI16**, C21orf63, RBP7, GABRP, CDC42BPG, SLC26A5, EFCAB4A, **FBLN5**, CBFA2T3, RGMA, **HOPX**, CASZ1, **ERBB3** |
| 上调（2） | **COL24A1**, **MMP1** |

生物学上连贯：**分化标志丢失（HOPX、ERBB3）+ 基质重塑（MMP1、COL24A1、FBLN5）**，与手稿主线 "stroma-dominated program" 自洽。

### 2.6 ✅ 已决策：趋势主口径 = 19（敏感性分析 = 1）

**用户 2026-10-18 拍板：采纳建议，正文主数字用 19，敏感性分析用 1。**

理由（复核确认）：模型 B（`gse78071_DE_ordinal_trend_pooledvar.csv`）与旧 pairwise 分析使用同一批 48 样本、同为 df 47.5、同一方差结构，唯一变量是把 `Mal−Ben` 换成线性 grade 斜率，因此 **62 → 19 是同口径可直接比较的真实差值**；模型 A（n=24, df 25.5）统计功效低、p 值偏保守，仅作敏感性分析。

**已执行的配套修改：**
- `scripts/python/08_figures.py` Fig 3B：原第三根柱读 `gse78071_DE_Trend.csv`（= Mal_vs_Ben 副本，62），现改为读 `gse78071_DE_ordinal_trend_pooledvar.csv`（= 19），x 轴标签改为 "Ordinal trend (continuous severity)"，并加注释说明旧表为代数恒等副本。
- Fig 4 改为从 `gse_pt_druggable_selected8_ruleA.csv` 读取规则选出的 8 靶点。

---

### 2.7 建议改稿文本（可直接用）

> **Ordinal grade analysis.** Because the `(borderline − benign) + (malignant − borderline)` contrast is algebraically identical to `malignant − benign`, we replaced it with a continuous severity regression (severity coded 1/2/3) fitted within the same linear model. **19 probes (19 genes) were associated with histological grade at FDR<0.05**, all nested within the 62 probes of the malignant-vs-benign contrast. The remaining 43 probes showed a binary malignant-vs-benign shift without monotone behaviour across the three grades. In a sensitivity analysis restricted to the phyllodes subset (n=24), one probe remained (ABO). Thus, at the transcriptome level, the dominant axis is the entity-level (PT vs fibroadenoma) contrast rather than a continuous grade-driven program.

**四条措辞纪律（重要）：**
1. 必须写明 **19 ⊂ 62**，否则读者会误算成 81 个基因。
2. **删除**"supporting the view that borderline is not an intermediate state" —— 线性趋势阴性**不排除**阈值型/阶梯型模式（Ben 低、(Bor,Mal) 同高）。要论证中间态需专门检验 `Bor − (Ben+Mal)/2 = 0`，当前未做。改为 "no evidence for a *linear* graded program"。
3. 必须给出**可检测下限**（中位方差探针 FDR<0.05 下可检出 ≈1.3 log2/grade 斜率），否则易被 "underpowered" 一句推翻。
4. 同一段不要混用 62/19/1 而不说明来源；固定为：正文 19，敏感性分析 1，并说明差异源于方差借用与自由度。

**ABO 处理：** 仅作脚注，不建生物学论点。理由：① 组织 ABO 转录本丰度受血型基因型强烈影响（O 型因移码显著偏低），恶性组仅 n=4，**单个 O 型个体即可驱动**；② GPL14951 上 ABO 仅 1 个探针，无内部重复；③ 它同时是 PT_vs_FA 显著基因（−0.87），更像追踪上皮/间质比例而非分级特异。三对比方向自洽（PT_vs_FA −0.870；Mal_vs_Ben −2.615；ordinal −1.225/grade；−2.615÷2=−1.307 与 −1.225 比值 1.07）可作为可信度脚注。

---

## 三、阻断项 B：可药靶 8 靶点表述（数字对、说法错）

### 3.1 手稿实际 8 靶点（此前内部沟通曾记错，以此为准）

**TGFB2、PDGFD、TWIST1、CCNE1、CDKN2A、COL1A1、MMP2、PARP1**
（CDK4/CDK6/FGFR1/PDGFRA **不在** 8 靶点内，它们只是可药靶表中 `druggable=True` 的普通行。）

### 3.2 核对结论

| 项目 | 手稿声称 | 实测 | 判定 |
|---|---|---|---|
| 候选基因总数 | 851 | 851 行（817 真基因 + 34 个 Agilent 探针 ID） | ⚠️ 应注明含未注释探针 |
| 优先靶点数 | "eight genes received at least one actionable agent" | 按此规则应得 **53**（`n_drugs>0`）；`druggable=True` 为 59 | ❌ 表述错误 |
| 8 靶点来源 | 隐含为筛选结果 | `rebuild_druggable.py:145` **手写硬编码** `MECH_VALID`(7) + PARP1 = 8 | ❌ 必须声明为人工提名 |
| DE 支持数 5/8 | TGFB2/PDGFD/TWIST1/CCNE1/CDKN2A | adj.P = 0.0174/0.0192/0.0251/0.0387/0.0385 | ✅ 一致 |
| 3 个不达标 | MMP2/COL1A1/PARP1 | 0.831/0.0541/0.2853 | ✅ 一致 |
| "five **independently** supported" | 独立证据层 | 候选池本身取自同一 DE 表，全表 **79.2%** 行 DE_sig=True | ❌ 循环论证，须删 "independently" |
| PARP1 无 PT 机制赋值 | 是 | CSV 中 PARP1 `pt_pathway=DNA_repair`、`pt_rel=True` | ❌ 与自身数据矛盾 |
| 字典条目 53 gene–agent–axis | 53 | 53 条目/53 基因，展开后 **111** gene–agent 对、81 唯一药名；其中 4 条无真实药物 | ⚠️ 口径应写清 |
| 8 靶点 Mal_vs_Ben 显著性 | 未说明 | **无一显著**（最小 adj.P=0.189, TWIST1） | ❌ 须明确 DE 支持来自 PT_vs_FA（实体对比）非分级 |
| 药名拼写 | `talazaparib` | 应为 **talazoparib** | ❌ 必改 |

### 3.3 需要你决策的一个科学问题

用可复现规则（`DE_sig(PT_vs_FA, adj.P<0.05) ∧ pt_rel ∧ n_drugs>0`）重筛，**同样得 8 个**，但成员不同：

| 入选（可复现规则） | adj.P | 字典药物 | 药物获批状态 | 手稿是否收录 |
|---|---|---|---|---|
| TGFB2 | 0.0174 | galunisertib | 未获批 | ✅ |
| PDGFD | 0.0192 | "PDGFR inhibitors(context)" | **无具名药物** | ✅ |
| TWIST1 | 0.0251 | "investigational(EMT)" | **无具名药物** | ✅ |
| CDKN2A | 0.0385 | "CDK4/6 inhibitors(context)" | **无具名药物** | ✅ |
| CCNE1 | 0.0387 | "CDK4/6 inhibitors(investigational)" | **无具名药物** | ✅ |
| **KIT** | 0.0421 | imatinib / sunitinib / avapritinib | **3 个均已获批** | ❌ 漏 |
| **TGFBR2** | 0.0429 | galunisertib | 未获批 | ❌ 漏 |
| **PDGFA** | 0.0493 | imatinib / nintedanib | **均已获批** | ❌ 漏 |

手稿用 COL1A1(0.054)、MMP2(0.831)、PARP1(0.285) **替换**了 KIT、TGFBR2、PDGFA。

### 3.4 ✅ 已决策：采纳方案 (a)（可复现规则）

**用户 2026-10-18 拍板：采纳方案 (a)。** 已实现并产出新表。

**执行结果**（`scripts/python/09_select_druggable_ruleA.py`）：

规则 `adjP_PT_vs_FA < 0.05 ∧ pt_rel ∧ n_drugs > 0` 从 851 个候选中**恰好选出 8 个**：

| 排名 | 基因 | 机制轴 | log2FC | adj.P | 药物 | 药物状态 |
|---|---|---|---|---|---|---|
| 1 | **TGFB2** | TGFb_ECM | +1.592 | 0.0174 | galunisertib; vactosertib | investigational |
| 2 | **PDGFD** | Stromal_other | −0.924 | 0.0192 | "PDGFR inhibitors(context)" | **无具名药物** |
| 3 | **TWIST1** | Stromal_other | +0.795 | 0.0251 | "investigational(EMT)" | **无具名药物** |
| 4 | **CDKN2A** | CellCycle_CDK | +1.279 | 0.0385 | "CDK4/6 inhibitors(context)" | **无具名药物** |
| 5 | **CCNE1** | CellCycle_CDK | +1.051 | 0.0387 | "CDK4/6 inhibitors(investigational)" | **无具名药物** |
| 6 | **KIT** | RTK_PDGFR_FGFR | −0.876 | 0.0421 | imatinib; sunitinib; avapritinib | **3 个均获批** |
| 7 | **TGFBR2** | TGFb_ECM | −0.719 | 0.0429 | galunisertib | investigational |
| 8 | **PDGFA** | RTK_PDGFR_FGFR | −1.049 | 0.0493 | imatinib; nintedanib | **均获批** |

**相比旧稿的改进：**
- DE 支持从 **5/8 → 8/8** 全部 FDR<0.05（规则本身即要求 DE 显著）
- 筛选**完全可复现**，审稿人跑脚本得同一结果
- 引入 **KIT、PDGFA** 两个带真实已批准药的靶点（旧稿的 PARP1 是唯一带已批药者，却 DE 最不显著）
- 主线叙事保住并强化：TGF-β 轴（TGFB2/TGFBR2）+ RTK 轴（KIT/PDGFA/PDGFD）+ EMT（TWIST1）+ 细胞周期（CCNE1/CDKN2A）

**⚠️ 必须如实写入正文的两条限制（不能回避）：**
1. **8/8 的 DE 支持全部来自 PT vs FA（实体对比）；在 Mal_vs_Ben 分级对比中 0/8 显著。** 不得表述为"随分级上调"或进展信号。
2. **仅 2/8（KIT、PDGFA）有已批准的直接靶向药**；4/8（PDGFD、TWIST1、CDKN2A、CCNE1）在 curated 字典中连具名化合物都没有，只是类别描述。因此 **"actionable agent" 一词不能用于全部 8 个**，建议改为 "eight targets meeting a reproducible selection rule (grade-independent differential expression + PT-mechanism annotation + curated candidate agent)"，并单独点明只有 KIT/PDGFA 具备已批准直接抑制剂。

**产出文件：**
- `gse_pt_druggable_selected8_ruleA.csv`（含 `agent_status` / `n_agents_approved` / `n_agents_investigational` 三列）
- `_druggable_ruleA_summary.txt`
- 脚本 `scripts/python/09_select_druggable_ruleA.py`

**已修正的字典错误：** `talazaparib` → **talazoparib**，已在 `scripts/R/07_annotate_drugs.R:58` 与 `scripts/python/rebuild_druggable.py:122` 两处源头改正。



---

## 四、通过项（可放心引用）

| 指标 | 手稿值 | 独立复算 | 判定 |
|---|---|---|---|
| GSE78071 n | 48 | 48（FA24 / Ben16 / Bor4 / Mal4） | ✅ |
| GSE231574 n | 33 | 33（Mal18 / Bor9 / Ben2 / FA2 / Meta2） | ✅ |
| GSE179458 n | 96 | 96（PT 子集 33） | ✅ |
| GSE235267 n | 16 | 16（FA4 / Ben6 / Mal6） | ✅ |
| PT_vs_FA DE | 737 | 737 探针（= 674 唯一基因；名义富集 2.66×） | ✅ 建议正文注明探针口径 |
| Mal_vs_Ben DE | 62（7 上调） | 62 探针 / 7 上调 / 55 下调 | ✅ |
| EPIC DMP（四对比） | 0/0/1/0 | 0/0/1/0（min adj.P 0.999999 等） | ✅ |
| DMR（DMRcate） | 0 | 表仅 1 行表头（0 数据行） | ✅ |
| 450k 探针数 | 361,964 | 361,964 | ✅ |
| 450k 显著 DMP | 0 | 0（min adj.P 0.1399） | ✅ |
| concordance overlap | 13,389 | 13,389（= 表内 18,544 symbol ∩ 双队列 DMP 基因） | ✅ 口径正确 |
| concordance observed | 51.2% | 0.5124 | ✅ |
| permutation null | 50.1% | 0.5013（95% CI 0.4929–0.5097） | ✅ |
| permutation p | ≈0.008（p<0.01） | 0.0045（种子 20260730，5 种子区间 0.0015–0.0080） | ✅ |

**说明：** 甲基化层的零结果手稿**原本就如实报告**（§3.1 明确写 0/0/1/0 DMP、450k 0/361,964），并配有显式效能分析，不存在过度陈述。

---

## 五、工程侧修复（已完成，非稿件问题）

| 项 | 状态 |
|---|---|
| `.env` DATASET_PATH 指向不存在的 `E:\CMCC Sync\...` | ✅ 已修正为 `E:\Medical Data\Dermatology\Phyllodes tumor` |
| `_env.py` 路径错时**静默回退到 Git 仓库** | ✅ 改为抛 RuntimeError（含修复提示），负向测试通过 |
| `_env.R` 路径错时**凭空造假目录树** | ✅ 加 `dir.exists` 断言 + stop()，验证未创建假目录 |
| 2 个历史补丁脚本内嵌已废弃错误值（`50.6%`/`p≈0.06`） | ✅ 归档至 `scripts/python/_archive/` 并写禁跑说明 |
| `08_figures.py` 硬编码 `C:/Users/tidus/...` | ✅ 改为基于 `__file__` 推导 |
| **8 个 R 脚本跨行 `else` 导致 parse 失败** | ✅ 已修（我的疏漏：只测了 `_env.R` 未测业务脚本），9 个脚本 `parse()` 全 OK |
| `_ordinal_trend.R` 误建在数据盘 | ✅ 已移入项目 `scripts/R/` 并改走 `_env.R`；数据盘 `scripts/` 已清除 |
| Biomni 执行后端（editable 指向已删路径） | ✅ 已修，不依赖 cwd 可 import，工具真实执行验证通过 |
| `current_data_version` = 2026-08-12 与实际不符 | ✅ 已改为 **2026-07-30** |

---

## 六、投稿前仍需补的占位符

| 位置 | 内容 | 责任方 |
|---|---|---|
| §8 Data Availability | 仓库 URL（项目目前**不是 git 仓库**） | 用户 |
| §9 Author Contributions | 13 个 CRediT 角色真名 | 用户 |
| 标题页 | 作者名 / 单位 / 通讯邮箱 / ORCID | 用户 |

---

## 七、开写建议顺序

1. **先定两件事**：趋势主口径（建议 19 主 + 1 敏感性）、可药靶口径（建议 (b) 人工提名）。
2. **可直接撰写**：Introduction、Discussion 论证、队列与样本量、PT_vs_FA 与 Mal_vs_Ben 结果、甲基化零结果与效能、concordance。
3. **需按本报告改写的两节**：§3 趋势段、§3.4 可药靶段。
4. **配图**：Fig 3B 的 "trend" 柱需替换为 19（或移除该柱）；Fig 4 无需改（8 靶点未变）。

---

*本报告由 Biomni 生物医学专家团复核产出，所有数字均经主理人独立复算，未修改任何数据文件。*
