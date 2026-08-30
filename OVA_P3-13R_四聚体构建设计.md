# OVA P3-13R 四聚体/多聚体可执行构建设计

## 1. 交付结论

本次没有在缺少 P3-13R 二聚体结构的情况下臆造新的 OVA 本体界面突变，而是把原有二聚体方案转换为可直接进入克隆评估的模块化候选：

1. **首选：T4-01-P3-pLI**
   
   ```text
   OVA P3-13R—(GGGGS)4—GCN4-pLI
   ```
   
   每条链含一个 OVA P3-13R；C 端 GCN4-pLI 预期驱动 4 条同源链形成平行四螺旋束，因此是本项目的四链四聚体候选。

2. **关键几何对照/备选：T4-02-WT-pLI**
   
   ```text
   OVA WT—(GGGGS)4—GCN4-pLI
   ```
   
   去掉 P3-13R 的已有自缔合突变，仅由 pLI 定义四聚态。如果 T4-01 出现更高阶聚集而 T4-02 主要为四聚体，说明 P3 界面可能在跨四聚体连接。

3. **探索性多价构建：M4-01-Tandem-P3**
   
   ```text
   OVA P3-13R—(GGGGS)4—OVA P3-13R
   ```
   
   如果两条串联链的两个 P3 模块都在同一对链之间闭合配对，可得到 4 个 OVA 模块；但它是“2 条链、4 个 OVA 模块”，不是四链同源四聚体，而且存在形成线性/网状多聚体的风险，因此不列为首选。

## 2. 输入序列复核

- 原 Excel 文件名写有“`OVA-P3-1P`”，但工作表中的实际标签与突变内容均为 **OVA P3-13R**；本次设计以工作表实际内容为准。
- OVA WT 与 OVA P3-13R 均为 **1158 bp / 386 aa**，源序列末端没有 stop codon。
- P3-13R 相对 WT 只有预期的 7 个氨基酸改变：`D13R; S68A; S76A; N78A; D190A; E191F; D192A`。
- 如果从起始 `ATG` 编码的 Met 开始计数，上述位置对应 `D14R; S69A; S77A; N79A; D191A; E192F; D193A`。
- 新候选中的 OVA 部分逐碱基保留源 Excel 序列；只在模块连接处加入 linker、四聚化域和终止密码子。

## 3. 为什么选择 GCN4-pLI

GCN4-pLI 的实验结构为平行、C4 对称的同源四聚体，PDB 编号为 **1GCL**；其 33 aa 序列为：

```text
RMKQIEDKLEEILSKLYHIENELARIKKLLGER
```

GCN4-pLI 的四聚态有晶体结构和经典实验依据，并曾作为融合蛋白的四聚化模块使用。因此它适合作为第一轮 proof-of-concept，而不是把“增加疏水性”当作四聚体设计。[RCSB PDB 1GCL](https://www.rcsb.org/structure/1GCL)、[Harbury 等，Science](https://doi.org/10.1126/science.8248779)、[GCN4-pLI 融合应用实例](https://pmc.ncbi.nlm.nih.gov/articles/PMC5391454/)

新增 pLI DNA 为通用哺乳动物偏好反向翻译：

```text
CGCATGAAGCAGATCGAGGACAAGCTGGAGGAGATCCTGAGCAAGCTGTACCACATCGAGAACGAGCTGGCCCGCATCAAGAAGCTGCTGGGCGAGCGC
```

## 4. Linker 设计

四聚化模块前使用 20 aa 柔性 linker：

```text
(GGGGS)4 = GGGGSGGGGSGGGGSGGGGS
```

DNA 使用不同的同义 Gly/Ser 密码子，减少精确核苷酸重复：

```text
GGTGGCGGAGGGTCCGGCGGTGGGGGAAGCGGAGGCGGTGGGTCTGGGGGTGGCGGATCA
```

该 linker 用于降低 OVA 球状结构与平行四螺旋束之间的空间冲突，但 linker 长度仍需通过结构预测和实验比较确认。

## 5. 构建规格

| 构建 | 架构 | AA/链 | CDS（含 TAA） | 单链理论 MW | 预期装配 | 预期装配 MW |
|---|---|---:|---:|---:|---|---:|
| D0-P3-13R | OVA P3-13R | 386 aa | 1161 bp | 42.78 kDa | 声称二聚体 | 85.55 kDa |
| D0-WT | OVA WT | 386 aa | 1161 bp | 42.88 kDa | 主要按单体解释 | 42.88 kDa |
| **T4-01-P3-pLI** | P3-13R—L20—pLI | 439 aa | 1320 bp | 48.06 kDa | 4 条同源链 | 192.23 kDa |
| **T4-02-WT-pLI** | WT—L20—pLI | 439 aa | 1320 bp | 48.16 kDa | 4 条同源链 | 192.65 kDa |
| M4-01-Tandem-P3 | P3-13R—L20—P3-13R | 792 aa | 2379 bp | 86.80 kDa | 2 条链/4 个 OVA 模块（理想情况） | 173.60 kDa |

理论分子量按未修饰多肽计算，不包含糖基化、磷酸化或其他翻译后修饰；实验表观分子量可与表中数值不同。

## 6. 首轮建议

首轮至少同时比较以下 4 个构建：

```text
D0-WT
D0-P3-13R
T4-01-P3-pLI
T4-02-WT-pLI
```

建议优先用 **SEC-MALS** 判断绝对分子量和均一性，并同时查看 SEC 峰形/回收率。判读重点：

- T4-01 是否主要集中在约 192 kDa 附近；
- T4-01 是否出现显著高分子量拖尾或 void-volume 峰；
- T4-02 是否比 T4-01 更接近单一四聚体；
- D0-P3-13R 是否确实接近约 86 kDa，而不是浓度依赖混合物。

只有在四聚态、单分散性和表达量都可接受后，才进入目标溶瘤病毒体系的功能评价。

## 7. 必须保留的风险判断

- **T4-01 是四聚体候选，不是已验证四聚体。** pLI 预期聚集 4 条链，但 P3 的已有界面可能在不同 pLI 四聚体之间继续配对，导致 8 聚体或更高阶网络。
- **T4-02 不只是普通阴性对照。** 如果目标优先级是“定义明确的四聚体”而不是“必须保留 P3 的全部突变”，它可能反而是更干净的最终架构。
- **M4-01 不是四链四聚体。** 它适合回答多价展示问题，但两个二聚位点也可能驱动聚合。
- GCN4-pLI 是酵母来源短肽，需评估新增免疫原性以及对目标病毒载量、表达和定位的影响。
- 当前候选未新增 signal peptide、跨膜区或纯化标签；这些应根据最终载体和表达定位另行决定。
- 所有完整 CDS 已加入末端 `TAA`。如果还要连接 C 端标签或其他模块，克隆时应删除该 stop。
- T4-01/T4-02 相对原始无 stop 的 OVA CDS 增加 **162 bp**（linker 60 bp + pLI 99 bp + stop 3 bp）；M4-01 增加 **1221 bp**，病毒载量评估时需计入。

## 8. 输出文件

- `OVA_P3-13R_四聚体候选序列.xlsx`：保留原工作表，并新增“`四聚体候选`”和“`模块与QC`”两个 sheet。
- `OVA_P3-13R_四聚体候选_AA.fasta`：所有对照和候选的完整氨基酸序列。
- `OVA_P3-13R_四聚体候选_CDS.fasta`：所有完整 CDS，均含末端 `TAA`。
- `build_tetramer_design.py`：可重复生成并自动核验上述序列。

## 9. 当前推荐决策

优先推进 **T4-01-P3-pLI + T4-02-WT-pLI** 成对测试。不要只做 T4-01；缺少 T4-02 时，一旦观察到高分子量组分，将难以区分它来自 pLI 四聚化还是 P3 界面导致的跨四聚体聚合。M4-01 仅在不接受外源四聚化短肽、且可以容忍多聚风险时作为补充路线。
