# OVA WT 与 OVA P3-13R 序列信息汇总

## 1. 文件概况

当前分析基于上传文件：

**`OVA WT及OVA-P3-1P突变体序列.xlsx`**

文件中包含多种 OVA（Ovalbumin，鸡卵清蛋白）相关 DNA 序列，其中重点核对的是：

- OVA WT
- OVA P3-13R

文件中 `WT 和 13R` 工作表明确标注了：

- `OVA WT DNA sequence`
- `OVA P3-13R DNA sequence`
- `Amino Acid Mutation`

P3-13R 对应的氨基酸突变为：

`D13R; S68A; S76A; N78A; D190A; E191F; D192A`

---

## 2. 这是什么序列？

这是 **DNA 编码序列（coding sequence / CDS / ORF）**，不是直接写出的氨基酸序列。

判断依据包括：

- 序列仅由 A、T、G、C 构成；
- 从起始密码子 `ATG` 开始；
- 可以按每 3 个核苷酸翻译成一个氨基酸；
- OVA WT DNA 长度为 **1158 bp**；
- 1158 = 386 × 3，因此恰好编码 **386 个氨基酸**。

---

## 3. OVA WT 是否为标准鸡卵清蛋白？

### 已确认

OVA WT DNA 翻译后得到 **386 aa** 的蛋白序列。

该蛋白序列与标准鸡卵清蛋白：

- 名称：Ovalbumin / OVA
- 物种：*Gallus gallus*
- UniProt：**P01012**
- 基因：**SERPINB14**
- 长度：**386 aa**

一致。

因此，文件中的 OVA WT 可以视为标准鸡卵清蛋白 OVA 的编码序列骨架。

### 名称注意

这里的 OVA 是：

**Ovalbumin（鸡卵清蛋白，鸡蛋清主要蛋白之一）**

不是：

**Chicken serum albumin（鸡血清白蛋白）**

二者不是同一个蛋白。

---

## 4. OVA WT 序列末端的特点

文件中的 OVA WT 长度为：

**1158 bp**

即：

`386 aa × 3 nt/aa = 1158 nt`

序列末尾没有额外的终止密码子，因此这里保存的是：

**完整 OVA 蛋白开放阅读框，但去掉 stop codon 的版本。**

如果加入一个终止密码子，例如：

`TAA`

则完整 CDS 长度会成为：

**1161 bp**

去掉 stop codon 常见于后续需要与其他蛋白、标签、linker 或其他功能模块进行融合表达的构建设计。

> 但仅凭“没有 stop codon”这一点，不能单独证明它一定用于融合蛋白。

---

## 5. OVA P3-13R 的氨基酸突变

文件给出的突变为：

| OVA 常用编号 | WT | Mutant | 变化 |
|---|---|---|---|
| D13R | D | R | Asp → Arg |
| S68A | S | A | Ser → Ala |
| S76A | S | A | Ser → Ala |
| N78A | N | A | Asn → Ala |
| D190A | D | A | Asp → Ala |
| E191F | E | F | Glu → Phe |
| D192A | D | A | Asp → Ala |

合计：

**7 个氨基酸突变。**

---

## 6. 氨基酸编号为什么与 DNA codon 编号差 1？

这是当前分析中非常重要的一点。

OVA DNA 从：

`ATG`

开始，对应蛋白 N 端：

`M...`

也就是说，DNA 的第 1 个 codon 编码起始 Met。

但是文件中的 OVA 突变编号采用的是一种 **不计最前端 initiator Met 的 OVA 编号方式**。

因此：

- 文件中的 D13
- 实际位于从 ATG 开始计数的第 14 个 codon

同理：

| 文件中的 OVA 编号 | 从 ATG 开始的实际蛋白位置 |
|---|---:|
| 13 | 14 |
| 68 | 69 |
| 76 | 77 |
| 78 | 79 |
| 190 | 191 |
| 191 | 192 |
| 192 | 193 |

后续进行：

- DNA 构建
- 定点突变
- 蛋白结构定位
- PDB/AlphaFold residue mapping

时，需要明确使用哪一套编号。

---

## 7. 7 个突变对应的 DNA 密码子变化

逐 codon 对齐得到：

| OVA 常用编号 | 从 ATG 起算 codon | DNA 位置 | WT codon | P3-13R codon | 蛋白变化 |
|---|---:|---:|---|---|---|
| D13R | 14 | 40–42 | `GAT` | `AGA` | D → R |
| S68A | 69 | 205–207 | `AGT` | `GCT` | S → A |
| S76A | 77 | 229–231 | `TCT` | `GCT` | S → A |
| N78A | 79 | 235–237 | `AAC` | `GCC` | N → A |
| D190A | 191 | 571–573 | `GAT` | `GCT` | D → A |
| E191F | 192 | 574–576 | `GAA` | `TTC` | E → F |
| D192A | 193 | 577–579 | `GAC` | `GCC` | D → A |

这些突变一共改变约 **13 个核苷酸**。

目前的逐位比较没有发现为了保持相同氨基酸而额外引入的大量同义密码子优化，因此 P3-13R 更像是：

**以 WT OVA DNA 为骨架进行定点突变得到的工程化变体。**

---

## 8. 190–192 区域是一个非常显著的局部变化

WT 的 OVA 190–192：

`D-E-D`

即：

- Asp
- Glu
- Asp

三个位点整体偏酸性、带负电。

P3-13R 中变成：

`A-F-A`

即：

- Ala
- Phe
- Ala

DNA 层面为：

### WT

`GAT GAA GAC`

翻译为：

`D   E   D`

### P3-13R

`GCT TTC GCC`

翻译为：

`A   F   A`

因此：

**DED → AFA**

这是整个变体里非常显著的一组理化性质改变：

- 去除了 3 个连续酸性侧链中的大部分负电特征；
- 引入 Ala；
- 中央引入一个芳香疏水残基 Phe。

这很可能显著改变该区域的：

- 表面电荷；
- 局部疏水性；
- 蛋白-蛋白接触性质；
- 局部构象或稳定性。

但仅根据序列不能证明这些变化具体形成了哪一种二聚界面。

---

## 9. 其他突变的理化性质变化

### D13R

`Asp (-)` → `Arg (+)`

属于明显的电荷反转。

可能显著改变该位置的：

- 局部静电环境；
- 盐桥可能性；
- 蛋白表面电荷分布。

### S68A、S76A

`Ser` → `Ala`

主要表现为：

- 去除羟基；
- 降低局部极性；
- 减少潜在氢键。

### N78A

`Asn` → `Ala`

主要表现为：

- 去除酰胺侧链；
- 降低局部极性；
- 改变氢键能力。

---

## 10. 关于“P3-13R 可以形成二聚体”

### 当前已知

提供该序列的信息来源声称：

**该工程化 OVA 可以形成二聚体。**

文件中的 7 个突变也确实会明显改变多个蛋白表面的：

- 电荷；
- 极性；
- 疏水性质。

其中尤其值得关注：

- D13R
- D190A
- E191F
- D192A

### 当前尚不能仅凭序列确认的内容

仅凭 Excel 中的 DNA 序列，目前还不能证明：

1. P3-13R 在溶液中一定稳定形成二聚体；
2. 二聚体的具体结合界面在哪里；
3. 哪一个突变是二聚化的关键驱动位点；
4. 二聚体是对称还是非对称；
5. 二聚化亲和力或 Kd 是多少；
6. 是否存在更高阶寡聚体；
7. 二聚化是否依赖浓度、pH 或盐浓度；
8. DED → AFA 是否直接位于二聚界面。

这些问题需要结合结构或实验数据判断。

---

## 11. 如何进一步验证二聚体

### 实验上

最直接的方法包括：

- SEC
- SEC-MALS
- AUC
- Native PAGE
- Mass photometry
- Analytical ultracentrifugation
- Cross-linking + MS
- Cryo-EM / X-ray（如果适用）

如果单体 OVA 约 43 kDa，则稳定二聚体的理论分子量约为：

**86 kDa**

SEC-MALS 能比较直接地区分：

- 单体；
- 二聚体；
- 更高阶聚集体。

### 结构计算上

可以考虑：

- AlphaFold 3 / multimer prediction
- AlphaFold-Multimer
- RosettaDock
- HADDOCK
- PISA interface analysis

重点需要检查：

1. 7 个突变是否聚集在同一个蛋白表面；
2. 两个 OVA 分子接触时这些位点是否位于界面；
3. WT 与 P3-13R 的界面能是否明显不同；
4. E191F 是否形成新的疏水 packing；
5. D13R 是否形成新的跨链盐桥；
6. 68/76/78 与 190–192 是否属于同一空间区域或协同界面。

---

## 12. 关于“用于溶瘤病毒”

### 当前已知

提供背景信息称：

**该改造 OVA 用于溶瘤病毒相关体系。**

因此可以确定这是一个与溶瘤病毒项目有关的 OVA 工程化序列。

### 目前不能仅凭该 Excel 判断

不能仅根据序列确认它在溶瘤病毒中的具体用途，例如不能直接确定它究竟是：

- 模型抗原；
- 免疫增强模块；
- 二聚化模块；
- 融合蛋白的一部分；
- 病毒表面展示模块；
- 分泌蛋白；
- 肿瘤免疫研究中的 OVA model antigen。

要确定这些，需要进一步获得：

- 病毒载体图谱；
- 表达 cassette；
- promoter；
- 上下游融合序列；
- linker；
- signal peptide；
- 文献或专利出处。

---

## 13. 关于 OVA 在肿瘤免疫研究中的背景

OVA 本身是免疫学中非常经典的模型抗原。

经典 OVA-derived T-cell epitope 包括：

**SIINFEKL（OVA257–264，常用 OVA 编号）**

它广泛应用于：

- OT-I T cell 模型；
- 肿瘤免疫研究；
- 抗原特异性 CD8+ T-cell response；
- 疫苗研究；
- 病毒载体抗原研究。

因此，如果一个溶瘤病毒项目使用 OVA，需要特别注意：

**OVA 可能既具有结构工程用途，也可能只是作为经典模型抗原使用。**

是否真的利用其“二聚化”作为主要功能，需要结合项目原始资料判断。

---

## 14. 当前最可靠的结论

截至目前，可以较有把握地总结为：

> 该 Excel 中的 OVA WT 是标准鸡卵清蛋白 Ovalbumin 的 DNA 编码序列，编码 386 aa OVA 蛋白，并去掉了终止密码子。OVA P3-13R 是在该 WT 骨架上引入 7 个定点氨基酸突变得到的工程化变体：D13R、S68A、S76A、N78A、D190A、E191F 和 D192A。对应 DNA 中约有 13 个核苷酸发生改变。最显著的局部改造之一是 OVA 190–192 的 DED → AFA，它会显著降低该区域的酸性并提高疏水性。背景信息称该变体具有二聚化性质并用于溶瘤病毒相关体系，但仅凭当前序列文件还不能确认具体二聚界面、二聚机制以及它在病毒载体中的具体功能。

---

## 15. 下一步建议

后续分析优先级可以是：

1. **找到 P3/P3-13R 的原始论文、专利或设计来源**
2. **把 7 个突变映射到 OVA 三维结构**
3. **预测 WT–WT 与 P3-13R–P3-13R 二聚体**
4. **计算 interface area / ΔG / buried SASA**
5. **确认 190–192 是否位于预测二聚界面**
6. **寻找溶瘤病毒载体图谱，确定 OVA 的具体表达形式**
7. **如果有实验数据，优先检查 SEC-MALS 或其他分子量证据**

---

# 16. 四聚体改造：当前设计理解

## 16.1 项目目标

当前新的工程目标是：

> 在现有 OVA P3-13R 二聚体基础上，进一步获得一个**稳定、定义明确、尽量均一的 OVA 四聚体**。

这里必须区分：

- **defined tetramer**：主要以四聚体存在，构型明确；
- **aggregation / higher-order oligomerization**：形成 4-mer、6-mer、8-mer 或更大的非均一聚集体。

真正希望获得的是前者，而不是简单通过增加疏水性让蛋白“聚得更多”。

## 16.2 当前最合理的总体思路：dimer of dimers

如果 P3-13R 已经通过一套既有界面稳定形成二聚体：

```text
A ── B
```

那么更自然的四聚体工程路线是：

1. 保留已有 P3-13R 二聚界面；
2. 在二聚体另一侧寻找新的可设计表面；
3. 设计一个独立的 **dimer–dimer interface**；
4. 让两个稳定二聚体进一步结合成四聚体。

概念结构：

```text
A ── B
│    │
│    │   ← 新设计的 dimer–dimer interface
│    │
C ── D
```

可以理解为：

```text
Interface 1 = P3-13R 已有 dimer interface
Interface 2 = 新设计的 dimer–dimer interface
```

理想情况下得到一个有限、封闭、对称的 **dimer-of-dimers homotetramer**，例如 D2 对称四聚体。

## 16.3 为什么不建议简单继续增加表面疏水性

P3-13R 已经存在显著的表面理化性质变化，例如：

```text
D190-E191-D192
        ↓
A190-F191-A192
```

即：

```text
DED → AFA
```

如果继续大量做：

```text
charged/polar → hydrophobic
```

虽然可能增强 self-association，但也容易产生：

```text
dimer
  ↓
tetramer
  ↓
hexamer
  ↓
octamer
  ↓
large aggregate
```

因此不能把“更容易聚集”等同于“形成稳定四聚体”。

四聚体设计的核心问题之一是：

> **怎样让组装在 4 个亚基处终止？**

所以应优先追求：

- 明确的几何互补；
- 有限的界面；
- 明确的对称性；
- 新旧界面的空间独立性；
- 避免可无限复制的 polymerization interface。

# 17. 四聚体设计的三条主要路线

## 17.1 路线 A：在 OVA 本体上设计第二界面

核心思路：

```text
P3-13R monomer
      ↓
existing dimer
      ↓
identify unused surface
      ↓
design dimer–dimer interface
      ↓
OVA tetramer
```

优点：

- 不需要引入较大的外源 tetramerization domain；
- 最终更接近“真正的工程化 OVA 四聚体”；
- 可以控制界面位置、对称性和构型。

缺点：

- 设计难度最高；
- 必须先知道 P3-13R 现有二聚体的真实或可信结构；
- 容易出现竞争构型；
- 需要严格的计算筛选与实验验证。

**这是长期目标最推荐的路线。**

## 17.2 路线 B：融合已知四聚化结构域

如果第一阶段核心问题是：

> “OVA 四聚化本身在目标溶瘤病毒体系里是否比二聚体更有价值？”

则可以先做：

```text
OVA-P3-13R
     │
   linker
     │
tetramerization domain
```

可考虑：

- 已知 tetramerization domain；
- tetrameric coiled-coil；
- 其他稳定形成 4-mer 的小型寡聚化模块。

优点：

- 实现四聚化相对直接；
- 适合做 proof-of-concept；
- 可先回答“四聚化有没有功能收益”。

缺点：

- 引入外源序列；
- 可能影响表达、稳定性、免疫原性和空间构象；
- 得到的是“OVA + tetramerization module”，不是 OVA 本体界面驱动的四聚体。

## 17.3 路线 C：二硫键或强相互作用锁定

可理论上考虑：

- 引入适当位置的 Cys；
- 形成 intermolecular disulfide bond；
- 或强化一个已经确定的四聚体界面。

但不建议在完全不知道四聚体构型时随机加 Cys，因为容易出现：

- 错误交联；
- 非均一多聚体；
- 表达/折叠问题；
- 更高阶聚集。

更合理的方式是：

> **先得到可信 tetramer geometry，再用二硫键作为后期锁定工具。**

# 18. 四聚体工程的关键前置问题

## 18.1 P3-13R 当前到底怎样形成二聚体？

目前知道 7 个突变：

```text
D13R
S68A
S76A
N78A
D190A
E191F
D192A
```

但仍不知道：

- 哪些位点直接位于 dimer interface；
- 哪些只是间接稳定突变；
- 二聚体的相对取向；
- 是否具有 C2 symmetry；
- P3 和 P3-13R 是否具有相同二聚构型；
- D13R 是否直接参与跨链盐桥；
- 190–192 的 AFA 是否直接参与 hydrophobic packing。

因此，四聚体设计前必须先建立可信的：

**P3-13R dimer structural model。**

## 18.2 7 个突变在三维空间上是否聚集？

线性序列中 13、68、76、78、190、191、192 相距较远，但折叠后可能：

1. 聚集在同一蛋白表面；
2. 共同组成连续的二聚界面；
3. 分成两块作用区域；
4. 部分位点只是通过稳定单体结构间接促进二聚。

这一点无法仅根据一维序列判断，因此应把 7 个突变映射到 OVA 三维结构上。

# 19. 推荐的计算设计流程

```text
OVA WT structure
        ↓
introduce P3-13R mutations
        ↓
P3-13R monomer validation
        ↓
P3-13R × 2
        ↓
dimer prediction / reconstruction
        ↓
map existing interface
        ↓
identify unused surface
        ↓
generate tetramer backbone / symmetry
        ↓
design interface sequence
        ↓
tetramer structure prediction
        ↓
interface scoring
        ↓
candidate ranking
        ↓
experimental validation
```

## 19.1 Step 1：确认 P3-13R 单体折叠

确认：

- 7 个突变不会明显破坏 OVA fold；
- 单体结构与 WT 基本一致；
- 特别关注 190–192 附近的局部结构变化。

## 19.2 Step 2：预测或重建 P3-13R 二聚体

可考虑：

- AlphaFold 3
- AlphaFold-Multimer
- RosettaDock
- HADDOCK

重点判断：

- 是否反复得到类似界面；
- 是否存在明显优势构型；
- 7 个突变是否位于预测界面；
- WT 是否也会得到同一界面；
- P3-13R 是否比 WT 更倾向该结构。

## 19.3 Step 3：现有二聚界面分析

建议分析：

- buried surface area / buried SASA；
- interface ΔG；
- shape complementarity；
- hydrogen bond；
- salt bridge；
- hydrophobic packing；
- buried unsatisfied polar atoms。

## 19.4 Step 4：寻找第二界面

理想的新界面应：

- 远离第一界面的关键残基；
- 不破坏 OVA core；
- 不产生过大的暴露疏水 patch；
- 能通过明确几何形成有限四聚体；
- 不产生无限 polymerization。

## 19.5 Step 5：构造 tetramer backbone

固定两个 P3-13R dimer，再搜索两个 dimer 之间的相对取向。

优先探索：

- D2 symmetry；
- C2 × C2 类型构型；
- 其他有限、封闭的四聚体几何。

## 19.6 Step 6：新界面序列设计

获得合理 backbone 后再设计 dimer–dimer interface。

可考虑：

- ProteinMPNN；
- Rosetta interface design。

设计原则：

- 尽量只修改新界面附近；
- 保护原有 P3-13R dimer interface；
- 保护 OVA core；
- 避免大面积暴露疏水 patch；
- 控制 mutation number。

## 19.7 Step 7：四聚体结构验证

候选序列重新进行：

- tetramer prediction；
- dimer prediction；
- monomer prediction。

目标：

- tetramer 构型稳定；
- dimer 中间体合理；
- monomer 仍保持 OVA fold；
- 不出现明显 misfold、alternative oligomer 或 polymeric assembly。

# 20. 四聚体候选筛选指标

不能只依据单一 AF confidence。

## 结构层面

- monomer fold preservation；
- tetramer symmetry；
- interface consistency；
- predicted confidence；
- alternative assembly tendency。

## 界面层面

- buried SASA；
- interface ΔG；
- shape complementarity；
- hydrophobic packing；
- hydrogen bonds；
- salt bridges；
- buried unsatisfied polar groups。

## 聚集风险

- exposed hydrophobic surface；
- self-association patch；
- polymerization geometry；
- higher-order assembly tendency。

## 工程可行性

- mutation number；
- expression；
- solubility；
- stability；
- fusion compatibility；
- virus expression constraints。

# 21. 四聚体实验验证

最终不能仅靠结构预测确认四聚体。

OVA 单体约 43 kDa，因此理论上：

```text
monomer  ≈ 43 kDa
dimer    ≈ 86 kDa
tetramer ≈ 172 kDa
```

优先建议 **SEC-MALS**，因为可以测绝对分子量，并区分：

- monomer；
- dimer；
- tetramer；
- higher-order aggregate。

其他可考虑：

- SEC；
- Native PAGE；
- AUC；
- Mass photometry；
- cross-linking；
- analytical ultracentrifugation；
- 条件允许时进行结构实验。

# 22. 项目阶段划分

## Phase 0：确认现有 P3-13R

- 重新确认序列；
- 确认表达；
- 确认真实 oligomeric state；
- 最好获得 SEC-MALS 数据。

## Phase 1：建立 P3-13R dimer model

- 找到可信的 dimer geometry；
- 确认 7 个突变的结构作用；
- 定义 interface 1。

## Phase 2：快速 tetramer proof-of-concept

可考虑：

```text
OVA-P3-13R–linker–tetramerization domain
```

目的是回答：

> 四聚化本身是否值得继续投入？

## Phase 3：intrinsic OVA tetramer design

如果四聚化确实有价值：

- 去掉外源 tetramerization module；
- 设计 OVA 自身第二界面；
- 获得 defined dimer-of-dimers tetramer。

## Phase 4：实验筛选

比较：

```text
WT
vs
P3
vs
P3-13R
vs
tetramer candidates
```

指标包括：

- expression；
- soluble fraction；
- SEC profile；
- SEC-MALS molecular weight；
- thermal stability；
- aggregation；
- 与目标溶瘤病毒体系相关的功能 readout。

# 23. 项目目录命名

推荐项目名：

```text
ova_p3_tetramer_design
```

推荐目录：

```text
ova_p3_tetramer_design/
├── 00_reference/
├── 01_p3_dimer/
├── 02_interface_mapping/
├── 03_tetramer_backbone/
├── 04_sequence_design/
├── 05_af3_validation/
├── 06_rosetta_analysis/
├── 07_candidates/
├── 08_experiment/
└── README.md
```

目录含义：

- `00_reference/`：WT、P3、P3-13R 序列，OVA reference structure，原始 Excel，文献/专利；
- `01_p3_dimer/`：P3/P3-13R dimer predictions、docking、AF3/AF-Multimer 输出；
- `02_interface_mapping/`：7 个突变位置、interface residues、SASA、电荷、结构可视化；
- `03_tetramer_backbone/`：dimer–dimer docking、symmetry generation、tetramer backbone candidates；
- `04_sequence_design/`：ProteinMPNN、RosettaDesign、mutation tables、designed sequences；
- `05_af3_validation/`：tetramer、dimer、monomer controls 和 confidence metrics；
- `06_rosetta_analysis/`：InterfaceAnalyzer、ΔG、buried SASA、shape complementarity、ranking；
- `07_candidates/`：最终 shortlisted candidates；
- `08_experiment/`：cloning、expression、SEC、SEC-MALS 和 functional assay data。

# 24. 当前项目的核心科学问题

目前可以概括为：

```text
标准鸡 OVA
    ↓
OVA P3 / P3-13R
    ↓
通过 7 个突变增强 self-association / dimerization
    ↓
确认真实 P3-13R dimer interface
    ↓
保留已有二聚界面
    ↓
设计独立 dimer–dimer interface
    ↓
defined OVA tetramer
    ↓
验证 tetramer 是否改善目标溶瘤病毒体系中的功能
```

当前最关键、仍缺失的信息是：

> **P3-13R 的真实二聚体结构和二聚界面。**

在这一信息确认之前，不建议直接大规模进行四聚体序列突变设计。

# 25. 当前最可靠的整体结论（更新版）

> 当前文件中的 OVA WT 是标准鸡卵清蛋白 Ovalbumin 的 DNA 编码序列，编码 386 aa OVA，并去掉终止密码子。OVA P3-13R 在该 WT 骨架上具有 7 个氨基酸突变：D13R、S68A、S76A、N78A、D190A、E191F 和 D192A；对应 13 个 DNA 碱基变化。其中 DED190–192 → AFA 是非常显著的局部电荷/疏水性质改变。背景信息称该变体可以形成二聚体并用于溶瘤病毒相关项目，但仅凭序列仍不能确认具体二聚界面和机制。
>
> 如果进一步将其工程化为四聚体，长期更合理的路线是首先明确 P3-13R 已有二聚界面，然后保留该界面，在另一侧设计空间上独立的 dimer–dimer interface，从而形成 defined dimer-of-dimers tetramer。短期也可以通过已知 tetramerization domain 做 proof-of-concept，以验证四聚化本身是否具有功能价值。所有四聚体设计都应优先避免形成无限聚合或非均一 aggregation，并最终通过 SEC-MALS 等实验确认真实 oligomeric state。
