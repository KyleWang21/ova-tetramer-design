# 中间文件清理候选（只读盘点）

盘点时间：2026-09-03 08:55 UTC；第一轮执行时间：2026-09-03 09:10 UTC。

## 建议优先清理（需先打包或确认）

1. `experiments/current_pending_v1_screen/ova_body_c4_06_ms25_r04/opendde/out/`
   - 约 992 MiB，10 个 OpenDDE 大型 `*_full_data_sample_0.json`。
   - 该候选已完成 OpenDDE 和 v1.0，但 required gate 失败；保留 `opendde_summary/`、`unified_final/`、日志和 10 个 CIF 后，原始大 JSON 可归档后删除。

2. 已完成 Accepted_AF3 批次中的 AF3 输入/置信度 JSON
   - 12 个已登记 AF3 批次合计约 865 MiB 的 `*_data.json`，另有约 0.6 MiB 的 summary JSON。
   - 23 条代表性 CIF、`full25_summary.tsv`/`seed1_strict_summary.tsv`、注册表和审计日志已在交付包或本地保留；这些 JSON 在不再需要重算时可压缩归档后删除。

## 可选清理

- 已完成批次中非代表模型的 CIF：约 1 GiB；当前交付包只需要每条候选 1 个排名 CIF。若仍需逐模型复核，则不要删。
- 已完成 Protenix/Rosetta 的冗余 PDB 和 decoy：Accepted_AF3 源目录约 215 MiB；应在确认 ProLIF、Rosetta 汇总和代表结构均已保存后再处理。

## 暂不清理

- E321 的 `run_next_*`、归档文件和日志：8 个任务仍在运行。
- E312 manifest/MSA、E319 AF3 输入目录：后续仍可能用于现有候选的 AF3 验证。
- `current_pending_v1_screen` 中除 `c06` 外的候选目录：大多仍处于 Protenix/化学计量/Rosetta/OpenDDE 后筛队列。
- `CURRENT_ALL_SAMPLES_STATUS.tsv`、Accepted_AF3 registry、各类 summary/TSV、研究日志、FASTA 和交付 CIF。

建议采用“先 tar.gz 归档、校验后再删除”的两步方式，不直接递归删除实验目录。

## 第一轮已执行

- 已将 166 个 Accepted_AF3 `*_data.json` 和 10 个 c06 OpenDDE `*_full_data_sample_0.json` 打包为 [`intermediate_json_archive_20260903.tgz`](deliverables/intermediate_json_archive_20260903.tgz)。
- 归档包含 176 个文件，原始总大小约 1.80 GiB；SHA-256：`299977983501ed632c523182e04ca95e6c8a96dd1be58852c2f69ddb9db01a70`。
- 已校验归档条目数、registry 中 23 个代表 CIF 路径和 c06 汇总文件；随后删除上述 176 个原始 JSON。
- 所有 CIF、汇总 TSV、registry、日志和正在运行/待验证目录均未删除；归档保留，可按需恢复。
