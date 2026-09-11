# 当前 OVA 四聚体样本与进度总览

更新时间：2026-09-03 14:43 UTC

详细逐样本总表见 [`CURRENT_ALL_SAMPLES_STATUS.tsv`](CURRENT_ALL_SAMPLES_STATUS.tsv)。该 TSV 是当前主线的单一总表，每行对应一个唯一的 386 aa OVA 序列；如果同一序列同时是 Accepted_AF3、E318 生成候选或 E321 锚点，相关身份合并在 `aliases`、`role` 和 `source_batch` 列中。

## 总计

| 类别 | 唯一序列数 | 当前阶段 | AF3 进度 |
|---|---:|---|---|
| Accepted_AF3 | 23 | 已完成严格五 seed × 五 sample | 575/575 模型完成；最高均值 ipTM 0.8660 |
| E312 冻结候选 | 11 | 等待火山 seed1 | 0/55 模型；每条计划 1 seed × 5 samples |
| E318 有效生成池 | 47 | 8/8 shard 已完成 | 0 个 AF3 模型；其中 14 条进入距离约束 shortlist，优先测试 |
| 历史纯 OVA 对照 | 4 | 旧批次结果 | 作为阴性/机制对照，不计入当前 Accepted_AF3 |
| **当前逐序列总表** | **85** | — | — |

E321 的 8 个锚点均来自上述 E318 或 Accepted_AF3 序列，因此没有新增唯一序列；其锚点身份在总表中合并记录。

## Accepted_AF3（23 条）

这些序列均已完成无模板 AF3 的 5 个独立 seed、每个 seed 5 个 samples，共 25 个模型，并通过当前严格的四链连通、零跨链原子 clash、双设计界面和二硫键约束。当前最高均值为 `ova_body_c4_26_ms25_r02` 的 0.8660（16 个突变，表面突变比例 0.8125）。

当前 `Accepted_v1.1=0`：

- 21 条仍在 ProLIF/Protenix-v2/化学计量/Rosetta 后筛队列中；OpenDDE不再阻塞这些候选；
- 2 条已完成 v1.1 核心后筛，但未通过全部 required gates：`ova_body_c4_073_ms`、`ova_body_c4_06_ms25_r04`；
- 没有任何一条可以标记为最终 v1.1 通过。

已经完整跑完后筛链的历史/当前样本逐阶段结果见 [`COMPLETED_POSTSCREEN_RESULTS.md`](COMPLETED_POSTSCREEN_RESULTS.md) 和 [`COMPLETED_POSTSCREEN_RESULTS.tsv`](COMPLETED_POSTSCREEN_RESULTS.tsv)。

## 实时资源与任务

- 中卫节点：E318 八个任务均已结束并完成归档；E321 八个任务均已获得资源运行，设计卡仍保持满载。
- 火山 `c20250508`：56/56 张 A100 已分配，当前空闲 0 张；E312 和其他 AF3 提交器严格等待空闲数超过 10。
- 本快照生成时仍没有 E318/E321 的 AF3 模型输出；所有 AF3 进度均以文件计数和严格汇总为准。
- **生成暂停**：已暂停提交新的 tied-AF2/序列设计批次。当前已启动的 E321 八个任务不强行取消，完成后只继续现有序列的归档和 AF3/后筛验证。

23 条的完整 AF3 均值、最低 ipTM、突变数、界面一致性和后筛状态在 TSV 中逐行给出；排名 CIF 仍在 [`deliverables/current_accepted_af3_ranked_20260903.tgz`](deliverables/current_accepted_af3_ranked_20260903.tgz)。

## E312：11 条冻结、尚未 AF3

E312 是 E311 八状态 tied-AF2 设计后按历史 AF3 标签和 Hamming 距离冻结的 11 条候选。它们均满足相对 P3-13R 不超过 24 个突变、表面突变比例不低于 0.80、天然 Cys 和 SIINFEKL 保持不变。

当前状态：火山队列 `c20250508`（`q-20251224105744-qhrgl`）空闲 GPU 低于提交门槛 `free>10`，所以 11 条尚未运行 seed1。表中的 `predicted_af3_mean_iptm` 和 `predicted_zero_clash_fraction` 是代理预测，不是 AF3 结果。

## E318：当前生成池

- 中卫任务：`ky-20260903-2300`–`2307`，Slurm `154006`–`154013`；
- 8 个 shard 已完成，最后一个（154013）已结束；归档已回收，正在由后处理脚本冻结候选并准备 E319 AF3 seed1；
- 已回收 8 个 shard 的轨迹状态和终点表，经约束过滤得到 47 条有效唯一序列；
- 其中 14 条通过历史标签加权和最小 Hamming 距离 3 的 shortlist，原计划进入 E319 无模板 AF3 seed1；当前生成已暂停，待资源和后续安排确认后再验证；
- E318 的 tied-AF2 ipTM 仅用于生成排序，不算 AF3 验收。

## E321：下一轮接力

E321 使用 E318 部分池提前启动，以免前 7 个完成 shard 释放的 A800 空闲：

- `ky-20260903-2500`–`2507`（Slurm `154107`–`154114`）均已运行（约 logits 12/14；`2507`约 logits 7/14）；
- E321 完成后不再启动新的设计批次，只做既有序列的 AF3 验证；
- 首次提交的 `154099`–`154106` 因 manifest 路径错误失败，已保留日志并修复后重新提交；
- E321 完成后仍必须经过无模板 AF3 五 seed × 五 sample 以及 v1.1 全部 required 筛选；OpenDDE仅作可选诊断。

## 其他对照

总表还保留 4 条历史纯 OVA 对照：原始 P3-13R 四链阴性基线，以及 `OVA-BODY-C4-C3HOT-AI01`、`OVA-BODY-C4-SS3-AI03`、`OVA-BODY-C4-P208F-AI01` 三条旧口径 C4-forming/机制对照。外接四聚螺旋的 `AI-C4-06` 不属于纯 OVA 序列，单独作为工程阳性对照，不纳入 85 条序列。

## 判定边界

“已产生”不等于“已验证”：

- AF2/tied-AF2、ProteinMPNN 和代理预测只能说明生成或排序信号；
- `Accepted_AF3` 才表示完整 25 模型 AF3 结构门槛通过；
- `Accepted_v1.1` 还要求 ProLIF、Protenix-v2、化学计量、Rosetta 等 required 项全部通过；OpenDDE不属于required；
- 当前最高 AF3 平均 ipTM 为 0.8660，目标约 0.90 尚未达到。
