# OVA-NC-C4-D2-92-AI01：非共价 OVA 四聚体设计报告

## 结论

`OVA-NC-C4-D2-92-AI01` 是当前主候选：386 aa 纯 OVA 本体，无 linker、binder 或外接四聚螺旋；未新增 Cys，只含 P3-13R 原有的 Cys12/31/74/121/368/383。AF3 五个独立 seed、共 25 个 C4 模型全部形成高置信四链连通装配，且 25/25 无 clash、25/25 无亚基间二硫键。

这是 AF3 计算设计命中，不是溶液实验已经证明的单一四聚体。强制六链输入在 seed3 出现 5/15 个通过模型，因此当前最准确的表述是“AF3 相对偏好且最稳健地形成 C4”，不能表述为绝对排斥 C6。

## 序列

```fasta
>OVA-NC-C4-D2-92-AI01|386aa|noncovalent_C4|no_engineered_Cys
MGSIGAASMEFCFRVFKELKVHHANENIFYCPIAIMSALAMVYLGAKDSTRTQINKVVRFDKLPGFGDAIEAQCGTAVAVHSSLRDIIDQITWHNRVYSFSLASNLYAEAGKMMDWEYWQCVEMLYEGGLQFINFQTAADQARILINRFVYLVTWRIIRNVLQPSSVDSQTAMVLVNAIVFRGLWEKAFKENATQAMPFRVTEQESKPVQMMYQIGLFRVASMASEKMKILELPFASGTMSMLVLLPDEVSGLEQLESIINFEKLTEWTSSNVMEERLIKVYLPKMKMEEKYNLTSVLMAMGITDVFSSSANLSGISSAESLKISQAVHAAHAEINTAGRYLLGEDWEAIDAASVSEEFRADHPFLFCIKHIATNAVLFFGRCVSP
```

相对 P3-13R 共 42 个突变：

`L88I,N89D,K93W,P94H,D96R,R105N,E110A,R111G,Y112K,P113M,I114M,L115D,P116W,L119W,K123E,E124M,R127E,E131Q,P132F,E144I,S148R,W149F,E151Y,S152L,Q153V,N155W,G156R,K182R,A191E,F192N,K278L,R285K,E337T,E341Y,V342L,V343L,S345E,A346D,E347W,A348E,G349A,V350I`

SIINFEKL 保持不变。

## 设计机制

1. 原始 P3-13R 在深 OVA MSA 条件下未通过 AF3 同源二聚体门槛，因此没有把原始低置信接触面当成已验证界面。
2. 用 AF2/BindCraft 风格的可微分同序列双链设计建立第一非共价界面，得到 `OVA-NC-DIMER-AF2-04`；其 AF3 C2 为 5/5 通过，ipTM 中位 0.66。
3. 保持该 AF3 二聚体不动，拟合其 179.59° C2 轴并搜索严格 D2 的正交轴方位与轴向中心。`phi=92°、offset=32 Å` 找到零 clash 的“二聚体的二聚体”设计骨架，并暴露连续的 C 端第二表面。
4. 在四链同序列约束下用 tied-ProteinMPNN 设计第二界面，重点形成 341–350 附近的非共价网络；最终仍由无模板 AF3 多 seed 独立判定，而不是用生成器分数定稿。

因此四聚化序列先验来自 OVA 表面的两类非共价界面，不来自亚基间二硫键或外挂模块。需要注意，AF3 没有原样恢复 MPNN 输入的严格 D2 刚体排布：允许链置换后，25 个 AF3 模型相对设计骨架的四链 Cα RMSD 中位为 19.13 Å；而 AF3 模型彼此相对最佳模型的中位 RMSD 为 0.77 Å。即 AF3 给出了内部一致但发生重排的紧凑四链构型，不能把 `phi=92°` 设计骨架当作最终结构真值。

## AF3 C4 结果

| seed | 模型数 | ipTM 中位 | 最弱链接口中位 | 通过 | 零亚基间 S–S |
|---:|---:|---:|---:|---:|---:|
| 1 | 5 | 0.70 | 0.69 | 5/5 | 5/5 |
| 2 | 5 | 0.66 | 0.65 | 5/5 | 5/5 |
| 3 | 5 | 0.82 | 0.79 | 5/5 | 5/5 |
| 4 | 5 | 0.82 | 0.77 | 5/5 | 5/5 |
| 5 | 5 | 0.79 | 0.75 | 5/5 | 5/5 |
| **合并** | **25** | **0.79** | **0.75** | **25/25** | **25/25** |

合并界面 PAE 中位 3.945 Å，最低 OVA 链 pLDDT 中位 86.02。最佳模型为 seed3/sample1：ipTM 0.83、最弱链 0.80、界面 PAE 3.30 Å。

## 强制链数对照

| 链数 | seeds / 模型 | ipTM 中位 | 最弱链接口中位 | 通过模型 | 判断 |
|---:|---:|---:|---:|---:|---|
| C2 | 1–3 / 15 | 0.65 | 0.65 | 15/15 | 稳定二聚装配中间体 |
| C3 | 1–3 / 15 | 0.23 | 0.19 | 0/15 | 失败 |
| C4 | 1–5 / 25 | 0.79 | 0.75 | 25/25 | 最稳健、最高总体置信 |
| C5 | 1–3 / 15 | 0.30 | 0.14 | 0/15 | 失败 |
| C6 | 1–3 / 15 | 0.26 | 0.33 | 5/15 | 只在 seed3 出现；弱于同 seed C4 |

同一 seed3 下，C4 的 ipTM/最弱链接口中位为 0.82/0.79，C6 为 0.65/0.66；因此 AF3 支持相对 C4 偏好，但不支持“绝对只有四聚体”的结论。

## 单体折叠与实验结构

C1 的 15 个模型 pTM 中位 0.91、最低链 pLDDT 中位 89.96。与 1.95 Å 未切割鸡 OVA 晶体结构 1OVA chain A 对齐：385 个共同 Cα 的全链 RMSD 中位 1.78 Å，高置信残基 RMSD 中位 1.25 Å，90.65% 的 Cα 在 2 Å 内。

主要偏差位于 OVA 原有的 70–82 柔性环及设计的 C 端 341–350 表面。后者在 C4 中接触复现率高，但局部 pLDDT 约 59–71，是当前最主要的实验风险点。

## 计算结论边界与实验建议

- AF3 支持该序列形成稳定非共价 C4，并相对 C3/C5/C6 更稳健；AF3 不直接给出溶液自由能、浓度依赖或唯一计量数。
- 首轮实验应优先做 SEC-MALS，并用 AUC 或 native MS 复核分子量；这是区分 C4 与可能的 C2/C6 亚群的关键。
- 用还原/非还原 SDS-PAGE 排除实际表达条件下意外链间二硫键；计算模型中 25/25 为零链间 S–S。
- 若实验出现明显 C6，可围绕 341–350 做负设计或逐位回退，但需重新进行 C4/C6 多 seed AF3 竞争，不能只看单个四链模型。

## 火山任务命名

火山文档与提交脚本均已强制使用 `ky-YYYYMMDD-NNN`，日期取北京时间，编号限定 `001–999`。本轮 AI 强化与链数对照实际使用 `ky-20260830-098` 至 `ky-20260830-106`；实验名、seed、shard 和标签均未写入任务名。
