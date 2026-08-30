# OVA P3-13R 纯本体四链设计结果

## 结论

已获得一个不含 linker、binder 或外接四聚螺旋的 386 aa OVA 本体序列：`OVA-BODY-C4-C3HOT-AI01`。它在 AlphaFold 3 的 8 个独立 seed、40 个模型中全部形成高置信、无 clash、四链连通复合物；39/40 同时形成化学干净的完整 12 条设计链间二硫键。因此它满足“能够在 AF3 中稳定形成四链复合物”的计算筛选目标。

必须同时保留一个重要限制：该序列的 C2/C3 输入也得到高置信结构，所以它是 **C4-forming 候选**，不是已经证明的“四聚体特异序列”，更不是溶液实验已证实的四聚体。

## 输入序列和阴性基线

- 源 Excel 中 P3-13R 为 386 aa，相对 WT 只有 `D14R/S69A/S77A/N79A/D191A/E192F/D193A` 七处变化。
- 深 OVA MSA、无模板条件下，原始 P3-13R 同源二聚体 25 个模型的中位 ipTM 为 0.17、界面 PAE 为 22.61 Å，未通过 AF3 二聚体判据。
- 同条件原始 P3-13R 四链输入 15 个模型的 ipTM/最弱链接口中位为 0.18/0.15、界面 PAE 为 26.44 Å，也未通过。
- 原始 P3-13R 单链折叠与 1OVA 实验结构一致：50 条 AF3 链的全链 Cα RMSD 中位 1.468 Å，高置信区域 0.627 Å。因此二聚失败不能归因于序列读错或单体整体塌陷。

权威表格：

- 二聚体基线：`experiments/e01d_dimer_msa/system_scores.tsv`
- 原始 P3 四链基线：`experiments/e29_parent_p3_c4_msa/system_scores.tsv`
- 1OVA 对比：`experiments/e23_experimental_structure_compare/summary.json`

## 最终 C4-forming 首选

`OVA-BODY-C4-C3HOT-AI01` 相对 P3-13R 的 19 个本体突变为：

`Q90E/K93E/P94C/D96C/V97N/S99Q/S101L/S152A/N155C/K182R/K190C/T194A/Q195R/A196C/Q210L/Y213C/Q214T/S237D/R340A`

其中 `P94C–Y213C`、`D96C–K190C`、`N155C–A196C` 是三组定向链间二硫键位点；其余位置来自成功/失败 AF3 状态、C4/C3 接触差热点和 tied-ProteinMPNN 多状态设计。SIINFEKL 与 N292 邻域保持不变。

### 八 seed AF3 结果

| 指标 | 结果 |
|---|---:|
| 独立 seed / 模型 | 8 / 40 |
| ipTM 中位 | 0.79 |
| 最弱链接口中位 | 0.78 |
| 界面 PAE 中位 | 1.723 Å |
| 最低链 pLDDT 中位 | 89.16 |
| 高置信四链连通 | 40/40 |
| 无 clash | 40/40 |
| 化学干净完整 12 键网络 | 39/40 |
| 天然 C74–C121 保留 | 40/40 |
| 非预期 S–S | 0/40 |

统一通过门槛为 `ipTM≥0.50`、最弱链接口 `≥0.45`、四链连通且无 clash；二硫键化学完整性单独统计，不用 ipTM 代替。

权威结果：

- `experiments/e52_c3hot_ai01_reseed/combined_8seed/combined_summary.json`
- `experiments/e52_c3hot_ai01_reseed/combined_8seed/per_seed_summary.tsv`
- `experiments/e52_c3hot_ai01_reseed/combined_8seed/combined_model_scores.tsv`

### 链数竞争

| 强制 AF3 链数 | 中位 ipTM | 最弱链接口 / pTM | 解释 |
|---:|---:|---:|---|
| C1 | — | pTM 0.92 | 单体折叠稳定 |
| C2 | 0.74 | 0.74 | 通过，存在强二聚倾向 |
| C3 | 0.75 | 0.73 | 通过，仍可形成三聚体 |
| C4 | 0.79 | 0.78 | 八 seed 稳健、当前最高 |
| C5 | 0.36 | 0.37 | 未通过；仅 4/15 形成完整设计键网络 |
| C6 | 0.34 | 0.40 | 未通过；0/15 形成完整设计键网络 |

这组结果中 C4 是最高分链数，且 C5/C6 明显失败，支持“C4 构象可稳定形成、不会简单无限延伸成更高聚体”。但 C2/C3 与 C4 的分差较小，仍不支持“只有 C4 能形成”。统一表格和图位于 `experiments/e53_c3hot_ai01_stoichiometry_extra/stoichiometry_summary.tsv` 与 `stoichiometry_competition.png`。

## AI 与热点设计过程

- 根据原始低可信二聚模型、SASA、pLDDT 和 MSA 保守性得到第一热点 `A99/A155/A157/A182/A184/A334/A336/A338`。其约 0.82 接触复现率只表示低可信模型中反复出现的表面，不表示已经验证的二聚界面。
- 真实运行 BindCraft 的 AF2 反向传播和 ProteinMPNN，得到可结合热点的 70 aa binder；全长 P3+binder 异源复合物通过 AF3，但六种融合四聚拓扑均失败，因此 binder 没有被冒充成最终 OVA 四聚序列。
- 纯非共价 OVA ProteinMPNN 两轮各 8 条候选均在独立 AF3 中失败。双二硫键设计能形成 C4，但存在双二聚拓扑和 seed 敏感性。
- 以 AF3 成功 C4 为正态、实际双二聚/C2/C3 为负态，进一步定位 `P94/Y213` 第三锁和 C4/C3 差分热点。最终首选由真实 AF3 构象驱动的 ProteinMPNN 生成，并由独立 AF3 多 seed 决定是否保留。
- 外接四聚螺旋的 `AI-C4-06` 只保留为工程阳性对照，未作为 OVA 本体四聚结果。

## 与实验 OVA 结构对比

最终候选 15 个 C4 模型的 60 条链与 1OVA 比较：全链 Cα RMSD 中位 1.518 Å，高置信区域 0.785 Å，92.7% 的匹配 Cα 在 2 Å 内。说明设计没有明显破坏 OVA 主折叠，但该比较不能替代表达、二硫键配对和溶液化学计量实验。

结果目录：`experiments/e45_c3_contact_negative_design/experimental_compare/`

## 交付文件

- 首选 AA：`OVA-BODY-C4-C3HOT-AI01_AA.fasta`
- 首选 CDS：`OVA-BODY-C4-C3HOT-AI01_CDS.fasta`
- 全部候选 AA/CDS：`OVA_P3-13R_四聚体候选_AA.fasta`、`OVA_P3-13R_四聚体候选_CDS.fasta`
- Excel：`OVA_P3-13R_四聚体候选序列.xlsx`
- 最佳 AF3 结构：`structures/OVA-BODY-C4-C3HOT-AI01_AF3_best.cif`
- C4 结构图：`figures/ova_body_c3hot_ai01_tetramer.png`
- 突变图：`figures/ova_body_c3hot_ai01_mutation_map.png`
- 八 seed 图：`experiments/e52_c3hot_ai01_reseed/combined_8seed/eight_seed_robustness.png`
- C4/C3 热点图：`experiments/e45_c3_contact_negative_design/hotspots/c4_vs_c3_hotspots.png`
- 完整实验日志：`research_log.md`

## 实验解释边界

AF3 支持该序列形成稳定四链结构，但由于 C2/C3 竞争同样较高，湿实验优先使用 SEC-MALS、AUC 或 native MS 确认分子量分布，并用非还原质谱/肽图验证三组设计二硫键。若实验目标是单一四聚体峰而非“可形成四聚体”，仍需继续做更强的链数负设计。
