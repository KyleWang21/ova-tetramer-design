# 当前最终候选产生路径总览

更新日期：2026-09-03

## 结论先行

当前有 **2 条主生成路径**，以及 **1 条备用/旧链路**：

1. **中卫 tied-AF2 多状态设计路径**：从当前最强 AF3 四聚体构象提取多骨架、多界面状态，在 8×A800 上生成新序列。
2. **火山 AF3 标签回灌主动学习路径**：用 E303/E307 的真实 AF3 标签、Hamming 核和保守上位性模型重组新序列，在无模板 AF3 中闭环筛选。
3. **链标号映射后的 bottom-k 救援路径（备用）**：从 E283/E284 的多骨架 tied-AF2 轨迹提取底部稳健状态，再做无模板 AF3；它等待主线前20条复核完成，不与主线争抢优先级。

这三条路径最后都汇入同一条 **v1.1 终审与交付路径**。ProLIF、Protenix-v2、C1/C2/C3/C5/C6 化学计量和 Rosetta 是共同的硬筛证据链；OpenDDE仅保留为可选第三模型诊断，不再阻塞终审。

## 路径 R-A：中卫 tied-AF2 多状态设计（主线）

### 生成与验证链

`E294 确定性强终点/当前 AF3 起点 → E311 8×A800 tied-AF2 → E312 无模板 AF3 seed1 → E313 Top4 五 seed×五 sample → v1.1 核心后筛`

- 输入是当前 Accepted_AF3 家族和冠军 `ova_body_c4_22_ms25_r01` 的四条独立 AF3 骨架。
- 每次更新同时约束四套骨架的两个真实界面状态（8 个状态），使用一条四链同序列参数共享更新。
- 只允许纯 386 aa OVA；相对 P3-13R 突变数 `<25`，表面突变比例 `≥0.80`，固定天然 Cys、SIINFEKL 和糖基化邻域，禁止新增 Cys 和亚基间二硫键。
- AF2/tied 分数只用于生成和排序；最终必须由无模板 AF3 的 5 seeds×5 samples 重新判定。

### 当前状态

- E311 的 8 个中卫任务 `ky-20260903-2200`–`2207`（Slurm `153894`–`153901`）仍在 8×A800 上运行，最近日志已进入 tied-AF2 的 soft 迭代，尚未回收归档。
- E311 完成后由 `finalize_e311_generation.sh` 回收并冻结 E312；随后 E312 才进入火山无模板 AF3。当前中卫 8 卡被 E311 保留；OpenDDE 观察器已停用。

### 权威证据

- 当前逐序列总表（Accepted_AF3、E312、E318部分池、E321锚点及历史对照）：`CURRENT_ALL_SAMPLES_STATUS.tsv`；概览：`CURRENT_ALL_SAMPLES_STATUS.md`

- `experiments/e311_c4_e294_detbest_continuation_tied_design/state_manifest.tsv`
- `experiments/e311_c4_e294_detbest_continuation_tied_design/targets/`
- `af3_pipeline/finalize_e311_generation.sh`
- `af3_pipeline/run_e312_after_e295_volc.sh`

## 路径 R-B：AF3 标签回灌/上位性主动学习（主线）

### 生成与验证链

`E303/E307 真实 AF3 标签 → E314 Hamming+保守 epistasis 重组 → E314 seed1 → E315 Top4 → E316 rank5–12 → E317 rank13–20 → v1.1 核心后筛`

- E314 从已测 AF3 标签和 354/356 局部碰撞标签中选择未测组合；模型只作为候选排序器，不替代结构验证。
- 每条候选先做 seed1 五 sample 严格检查；通过者再做 seed2、seed3、seed4，三轮都通过才做 seed5。
- 这里的 `seed2` 指 AF3 的第二个独立随机种子；该种子下仍运行 5 个 diffusion samples，因此 `seed2` 是 1 个 seed / 5 个样本，不是单独的 1 个样本。完整 seed1–5 复核为 5×5=25 个模型。
- 硬门槛是 25/25 模型、平均 ipTM、最低 ipTM、四链 C4 连通、零跨链原子 clash、双设计界面、天然链内二硫键完整、零亚基间二硫键和单体结构保持。
- 火山只使用队列 `q-20251224105744-qhrgl`（`c20250508`），仅当空闲 GPU `>10` 时提交，单批最多 8 张 A100；任务名统一为 `ky-YYYYMMDD-编号`。

### 当前状态

- E314 seed1：32 条中 20 条严格通过。
- E315：`ova_body_c4_26_ms25_r02`（16 突变，平均 ipTM `0.8660`）和 `ova_body_c4_13_ms25_r03`（16 突变，平均 `0.8564`）进入 Accepted_AF3。
- E316：8 条中 4 条新增 Accepted_AF3：`c02`、`c07`、`c28`、`c21`。
- E317：排名 13–20 的 8 条已完成 `seed2→seed3/4→seed5` 的逐seed fail-closed 复核；5 条（`c30`、`c09`、`c16`、`c06`、`c29`）完成 25 模型并通过，3 条因碰撞或未满足逐seed门槛淘汰。当前 Accepted_AF3 总数为 **23**，其中最高完整25模型均值为 `c30=0.8624`（16突变，表面比例0.875）。E317权威汇总为 `experiments/e317_e314_rank13_20_af3_multiseed/full25/full25_summary.tsv`。
- E311：中卫 8×A800 的四骨架八状态 tied-AF2 续设计已 8/8 完成并回传；生成池 40 条（13–20 突变），经 Hamming 去重和历史标签重排后冻结 E312 的 11 条（13–19 突变，表面突变比例 0.875–1.0）。
- E312：已生成 `CANDIDATES_FROZEN`，将对这 11 条分别运行无模板 AF3 seed1×5 samples；当前 c20250508 占用 56/56 张 A100（free=0），因此尚未提交，watcher 将在 free>10 时自动提交。E312 不再等待旧 E295/E287 备用链。
- E318：在 E311 完成后立即接力占用中卫 8×A800；4 个 E311 新轨迹锚点与 4 个高均值 Accepted_AF3 锚点组成八状态 tied-AF2 设计，任务 `ky-20260903-2300`–`2307`（Slurm 154006–154013）已提交。前7个shard已完成，最后一个154013目前处于 hard 2/4；远端归档watcher已恢复，完整归档后由 `af3_pipeline/finalize_e318_generation.sh` 冻结新序列并接入 E319 无模板 AF3 seed1，随后由 E320 对严格命中者做 seed2–5 五 seed×五 sample 复核。E318 的 AF2/tied 分数不计入 Accepted_AF3。
- 为保持中卫持续产出，E321 已利用前7个完成shard的部分候选池提前接力：8个锚点已冻结，`ky-20260903-2500`–`2506`（Slurm 154107–154113）运行，`2507`（154114）等待154013释放的A800。E321结果仍须完整轨迹、无模板AF3和v1.1核心后筛；E322/E323自动衔接脚本继续等待E321归档。

### 权威证据

- `experiments/e314_c4_e307_epistasis_reflow_af3_seed1/seed1_strict_final/`
- `experiments/e315_e314_top4_af3_multiseed/full25/full25_summary.tsv`
- `experiments/e316_e314_rank05_12_af3_multiseed/full25/full25_summary.tsv`
- `experiments/e317_e314_rank13_20_af3_multiseed/selected_candidates.tsv`
- `af3_pipeline/run_e316_e317_remaining_e314_multiseed.sh`
- `experiments/current_accepted_af3_registry/accepted_af3_registry.tsv`
- `af3_pipeline/finalize_e318_generation.sh`（E318→E319→E320 自动衔接）
- `af3_pipeline/submit_next_zhongwei_tied.sh`、`af3_pipeline/finalize_e321_generation.sh`（持续生成接力）

## 路径 R-C：chainmapped bottom-k 救援（备用，不计入当前主线数量）

### 生成与验证链

`E283 接触图重映射 → E286 底部多状态 tied-AF2 → E287 无模板 AF3 seed1 → E288 Top4 多 seed`

- 该链路从经过接触图 Jaccard 映射审计的四骨架八状态中提取 bottom-k 稳健离散状态。
- 它保留了旧批次中尚未被主线标签回灌覆盖的序列多样性，作为主线失败或达到平台期后的补充来源。
- 当前 `run_e287_after_e284_volc.sh` 等待 E314 的 `TOP20_MULTISEED_VALIDATION_COMPLETE`，因此尚未提交新的火山任务。
- E295/E296 等确定性强终点旧链路也属于历史/备用验证，不作为当前两条主线之外的独立成功证据。

## 共用的最终验收与交付路径

任何生成路径的序列都必须经过同一条 fail-closed 终审：

`无模板 AF3 5×5 → ProLIF 界面化学 → Protenix-v2 多 seed → C1/C2/C3/C5/C6 化学计量 → Rosetta 200 次 FastRelax → v1.1 单一总表 → 排名 CIF/FASTA/总结`

当前 Accepted_AF3 只表示 AF3 25 模型结构门槛已通过，不表示已经通过全部 v1.1 后筛。`TOP20_MULTISEED_VALIDATION_COMPLETE` 已生成，23条候选已进入统一 v1.1 队列；ProLIF/Protenix/Rosetta按fail-closed required统计，OpenDDE只作为可选诊断，不能把 Accepted_AF3 直接称为 Accepted_v1.1。

## 不再计入最终路径的历史方法

- **BindCraft binder/OVA 融合**：真实运行过，但融合拓扑未稳定形成目标纯 OVA 四聚体；只保留为 AI 方法学对照。
- **GCN4-pLI/CC-Tet 等外接四聚螺旋**：可作为工程阳性对照，但不符合纯 OVA 本体目标。
- **引入亚基间二硫键的四聚体设计**：已被当前目标禁止，不进入最终候选。
- **单一固定 C4 骨架的 ProteinMPNN 分数排序**：曾用于生成，但无法独立证明无模板 AF3 四聚体；现在只作为序列兼容性或多骨架排序信号。

因此，当前报告中的“路径数量”应写成：**2 条主生成路径 + 1 条备用生成路径 + 1 条所有路径共用的终审/交付路径**。
