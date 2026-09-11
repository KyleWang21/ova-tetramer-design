# OVA-LM-C4-17-B 最终计算结果

## 结论

`OVA-LM-C4-17-B` 是当前首选候选：386 aa 纯 OVA，同源四链使用完全相同序列；相对 P3-13R 共 17 个突变，无外挂模块、无新增 Cys、无设计链间二硫键，SIINFEKL（258–265）保持不变。

AF3 无模板 seeds1–5、每个 5 samples，共 25/25 模型通过预设四聚体置信度门槛，且全部四链连通、无 clash、零链间 S–S。独立 Protenix-v2 seed15 也支持完整非共价四链网络。这里的“形成四聚体”仍是计算预测，最终需用 SEC-MALS、AUC 或原位质谱验证溶液计量数。

## 序列与突变

突变：`L88I,N89D,K182R,A191E,S271W,K278L,R285K,E337T,E341N,V342L,V343L,S345E,A346D,E347G,A348G,G349A,V350I`

```text
MGSIGAASMEFCFRVFKELKVHHANENIFYCPIAIMSALAMVYLGAKDSTRTQINKVVRFDKLPGFGDAIEAQCGTAVAVHSSLRDIIDQITKPNDVYSFSLASRLYAEERYPILPEYLQCVKELYRGGLEPINFQTAADQARELINSWVESQTNGIIRNVLQPSSVDSQTAMVLVNAIVFRGLWEKAFKEFATQAMPFRVTEQESKPVQMMYQIGLFRVASMASEKMKILELPFASGTMSMLVLLPDEVSGLEQLESIINFEKLTEWTSWNVMEERLIKVYLPKMKMEEKYNLTSVLMAMGITDVFSSSANLSGISSAESLKISQAVHAAHAEINTAGRNLLGEDGGAIDAASVSEEFRADHPFLFCIKHIATNAVLFFGRCVSP
```

天然 Cys 位置仅为 12/31/74/121/368/383。17 个突变中，182/191/271/278/337/341–350 共 14 个在 AF3 通过模型中反复位于链间界面；88/89/285 更偏外围。

## AF3 无模板多种子

| 指标 | 结果 |
|---|---:|
| seeds × samples | 5 × 5 = 25 |
| 通过置信门槛 | 25/25 |
| ipTM 中位数 | 0.74 |
| 最弱链接口中位数 | 0.71 |
| 接触链对 PAE 中位数 | 4.57 Å |
| 最低链平均 pLDDT 中位数 | 88.74 |
| 四链全连通 / 无 clash | 25/25 / 25/25 |
| 链间 S–S | 0/25 |
| 四条天然 C74–C121 几何完整 | 24/25 |

唯一化学检查例外是 seed4 sample1 未把其中一个天然链内 C74–C121 放到 2.6 Å 内；该模型仍无链间或非预期 S–S，序列也未改变天然 Cys。

## Protenix-v2 独立复核

设置：464M、无模板、1024-row 候选专属 MSA、10 cycles、200 diffusion steps、OpenFold LayerNorm、seed15。

结果：pLDDT 86.21、pTM 0.699、ipTM 0.615；六个链对的 8 Å 残基接触数为 16/43/43/39/47/17，四链全连通、无 clash；四条天然 C74–C121 均完整，零链间 S–S。该单 seed 结果是独立支持，不替代 AF3 多 seed 统计。

## 实际设计流程

1. 旧 42 突变 C4 的 AF3 构象提供两类链间界面的空间起点。
2. 完整 4×386 AF2 反向传播需约 279–300 GiB，单张 80 GB 卡不可运行，因此没有宣称完成四链联合反传。
3. 改用同一 tied 序列的两个全长二聚状态：AB 代表第一界面，AC 代表第二界面；每个状态同时约束单体折叠接触和目标链间接触，并按相对 P3-13R 的全序列突变数惩罚。
4. 冷启动与 warm-start 稀疏化并行；warm start 从已验证界面序列开始删除冗余突变。16 张 A100/A800 共生成 28 条不重复的 ≤25 突变候选，20 条进入 AF3，其中 15 条 ≤20。
5. 最终完整四链、无模板 AF3 是唯一入选标准。`17-B` 因 25/25 跨 seed 通过而优于总体分数略高、但 seed5 全部失败的 `17-A`。

AF2 在这里是可微分的界面稀疏器和候选生成器，不是四聚体成功证据；完整四聚状态由 AF3 多 seed 和 Protenix-v2 独立复核判断。

## 主要文件

- `OVA-LM-C4-17-B.fasta`：最终序列
- `OVA-LM-C4-17-B_best_AF3.cif`：最佳 AF3 模型（seed2 sample1）
- `combined_model_scores.tsv`：25 个 AF3 模型逐项评分
- `combined_summary.json`：AF3 汇总
- `per_seed_summary.tsv`：逐 seed 汇总
- `mutations.tsv`：突变与保护位点
- `interface_hotspots/interface_positions.tsv`：界面复现率
- `protenix/seed15_score.json`：独立 Protenix-v2 评分
