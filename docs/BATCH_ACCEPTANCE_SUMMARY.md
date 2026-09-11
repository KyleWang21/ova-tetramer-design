# OVA 非共价四聚体历次批次验收汇总

更新时间：2026-09-03 UTC

## 状态定义

- **生成结果**：AF2 可微分设计、ProteinMPNN 或代理模型产生的序列；尚无 AF3 验证。
- **严格 seed1**：一个无模板 AF3 seed 的 5 个 samples 全部通过冻结的坐标、拓扑、二硫键、双设计界面和单体折叠门槛。
- **Accepted_AF3**：五个独立 AF3 seeds、每个 5 samples，共 25 个模型通过完整 AF3 门槛。
- **Accepted_v1.1**：在 Accepted_AF3 基础上，再通过序列约束、ProLIF、Protenix-v2、化学计量和 Rosetta 等全部适用 required 门槛；OpenDDE仅作为可选诊断，不再参与淘汰。

单一高 ipTM、单一最佳模型、AF2 接口分数或严格 seed1 均不能替代后两级验收。

## 批次汇总

| 批次 | 设计规模/突变 | AF3 结果 | 其他模型与后筛 | 当前结论 |
|---|---:|---|---|---|
| E55–E103 早期非共价探索 | 多种 D2、AF2 梯度和 MPNN 路线 | 多数在无模板 AF3 seed1 或独立 seed 中失败 | 未进入完整后筛 | 仅保留失败轨迹和热点信息，不保留最终序列 |
| E104 `OVA-NC-C4-D2-92-AI01` | 42 突变 | 当时口径下 25/25 四链连通、零 clash；ipTM 中位 0.79 | Protenix-v2 三 seed ipTM 0.301–0.444，只支持弱四链网络 | 超出当前 ≤24 突变预算，且跨模型支持不足；降级为设计母本/对照 |
| E110–E114 `OVA-LM-C4-17-B` | 17 突变 | 旧口径 25/25 通过；ipTM 中位约 0.74 | 单 seed Protenix-v2 ipTM 0.615 | 证明低突变非共价 C4 可行，但未达到后续 ipTM>0.80 目标，不是 v1.0 Accepted |
| E117–E144 旧 Final20 | 20 条，20–21 突变为主 | 25 模型平均 ipTM 0.806–0.8468；按 v1.0 重算只有 4–22/25 模型零 clash，0/20 Accepted_AF3 | ProLIF 16/20 通过；化学计量 20/20 通过；Protenix 0/20、Rosetta 0/20；OpenDDE 无合格代表 | **0/20 Accepted_v1.0**。序列表面突变比例仅 0.684–0.750，20/20 也不满足当前 ≥0.80 门槛 |
| E158/E165 `OVA-C4-073` | 17 突变；表面 14/17=0.824 | **Accepted_AF3**：25模型平均ipTM 0.8504、最低0.84；25/25零clash、完整C4、天然链内二硫键、零链间二硫键、双设计界面 | Stoichiometry通过：C4 margin 0.3624、第二界面25/25相对C2独有；Protenix-v2仅1/5几何干净且均值0.6771，OpenDDE 0个可选代表；Rosetta 200/200保持C4和零clash但最弱H-bonds中位0、unsat 9.03/1000 Å²。 | sequence/AF3/ProLIF/stoichiometry通过，Protenix/Rosetta/OpenDDE失败；保留为设计母本，**不能成为Accepted_v1.0** |
| E159 `OVA-C4-105` | 17 突变 | 25 模型平均 ipTM 0.7128；seed2 均值 0.278；14/25 零 clash；Accepted_AF3=0 | 未进入正向完整后筛 | 随机 seed 敏感负标签 |
| E162 回退库 | 16 条、13–17 突变 | 10/16 严格 seed1；优先复核的 15、14、13、14 突变体均在独立 seed2 的零-clash门槛失败 | fail-fast，未进入后筛 | 仅保留为 AF3 标签；0 条新增 Accepted_AF3 |
| E168 `OVA-C4-073`低突变回退修复 | 8条、14–16突变 | 5/8严格seed1。REV002完整25模型平均/最低ipTM 0.8516/0.84但只有23/25零clash；REV010 seed2均值0.328；REV003通过seed2但seed3只有4/5零clash；REV006/REV001 seed2分别只有4/5、3/5零clash | 五条命中均按逐seedfail-closed淘汰 | **0条新增Accepted_AF3**；15–16突变回退尚未复现母本073的25/25零clash |
| E174 六亲本等位组合 | 16 条、17–20 突变 | 16/16 完成；其中 7/16 严格 seed1，命中者平均 ipTM 0.830–0.852 | 09/02/06独立seed2淘汰；01和03完成25模型但分别仅20/25零clash；15完整25模型平均ipTM 0.736且19/25零clash；11为唯一25/25严格通过者，平均ipTM 0.8424。11的AF3-medoid ProLIF通过；Protenix-v2五seed平均ipTM/pTM为0.6950/0.7579且仅1/5几何干净。Rosetta 200/200保持拓扑和零clash，但最弱界面氢键中位0、未配对极性原子9.77/1000 Å²。修正MSA后的OpenDDE十seed平均/最大ipTM仅0.2694/0.3323，0个几何可选代表 | **新增1条 Accepted_AF3：`AI-MULTIHIT-11`**；其余6条严格seed1命中均被多seed复核淘汰。11在Protenix、Rosetta和OpenDDE均失败，不能成为Accepted_v1.0 |
| E252/E253/E289/E290/E291 二阶等位AI与AF3碰撞修复 | E252 32条、16–21突变；E289修复32条、14–20突变；表面比例均≥0.80 | E252候选09在seed1–4共20模型平均/最低ipTM 0.8495/0.84，但seed3/4各仅4/5零clash。E289用实际碰撞位点、Accepted等位基因和双目标Hamming AI产生修复库，seed1严格命中9/32。E290的19突变`E348A`修复体通过25模型：平均/最低ipTM 0.850/0.83，25/25零clash、完整C4、天然链内二硫键、零链间二硫键、双设计界面 | 新候选ProLIF：12/19突变在任意界面、6/19非VdW，主热点四链对称；两界面非VdW突变数1/3，故ProLIF失败。其余v1.0阶段准备中。E291继续复核另外5条seed1命中 | **新增1条Accepted_AF3：`E289-E348A-06`/`ova_body_c4_06_ms25_r04`**；Accepted_v1.0暂为0 |
| E155–E276 新生成池 | 12–24突变，均按表面比例≥0.80过滤 | E177的2条和E178的04/15严格seed1命中均在独立seed淘汰；E212/E216及后续批次等待新AF3 | E219冻结E220的4条。E221跨AF3–Protenix tied-AF2八分片完成：8条原始有效序列均15–20突变、表面比例1.00，经去重/距离约束冻结E222共7条；与E223快速MPNN联合后冻结E224 Top32，为16–22突变、表面比例0.85–0.947。审计发现旧slim archive遗漏trajectory；从E217/E219/E221补回414条原先未进入冻结池的有效近离散序列，E257冻结Top32。E235八分片完整回收后，从113条有效唯一AF2状态和222条有效代理序列中冻结E236 Top32；另发现14突变、双状态AF2 ipTM 0.739/0.637的直接轨迹被代理排到192名，故E259独立冻结按实测最弱tied-AF2 ipTM排序的8条救援序列（12–16突变、表面比例1.00）。E246/E247完整MPNN合并1278条有效唯一序列并冻结E248 Top32（12–20突变、表面0.813–0.944）。E250二阶等位对ridge与Hamming核按40:60混合，冻结E252 Top32并以Top4启动E251。E243八任务完成并通过完整轨迹归档；207条有效生成短名单冻结E244 Top32，为13–20突变。E243最佳完全离散状态为19突变、表面1.00、双状态AF2 ipTM 0.687/0.610；已提交E261任务1152–1159。E251八任务完整回收后冻结E254 Top32；E268实际选择soft=1.00、18突变、表面1.00、双状态0.690/0.672起点，任务1224–1231已提交，结果冻结E269。E265则在E261四context结果中构建≤20突变跨骨架共识，结果冻结E266。E271进一步冻结512条真实序列，在4套AI11-AF3与4套OVA073-Rosetta四聚体骨架上用精确tied-ProteinMPNN条件NLL、每骨架4种解码顺序独立重评分；按最差/平均骨架NLL及跨家族稳定性与AF3历史标签联合筛选后冻结E272。E261出现`0.715/0.682`连续态后另预注册E274，保持E265 hard规则不变，独立将最佳soft≥0.65 argmax起点跨四骨架长退火稳定，结果冻结E275/E276 | E251完成的卡已由E261即时补位，中卫保持8/8；E261后自动衔接E268/E265/E271/E274。火山按`free>10`门槛等待；唯一串行顺序为E252→E257→E259→E236→E254→E244→E248→E262→E264计量补测→E266→E269→E272→E275→旧批次。严格seed1命中均自动补齐25模型。所有AF2/MPNN/代理绝对值均不替代AF3验收 |

## 当前总验收结论

- 当前注册表中的 **Accepted_AF3**：23 条；完整逐样本状态（含E312、E318、E321）见 `CURRENT_ALL_SAMPLES_STATUS.tsv`。
- 历史首批明确命名的 Accepted_AF3 母本包括 `OVA-C4-073`、`E289-E348A-06`（注册名`ova_body_c4_06_ms25_r04`）和 `AI-MULTIHIT-11`，现已与后续批次统一纳入23条注册表。
- 新增 **Accepted_v1.1**：0 条（当前 required 集合不含 OpenDDE）。
- 旧 Final20 经核心 required 阶段追溯重算：0/20 全阶段通过，不能继续称为最终验收成功。
- 当前 E174：7 条严格 seed1 命中中，09/02/01/03/06/15六条在独立多seed复核失败；只有11通过25模型，新增1条多seed通过者。
- 目前最高可信五 seed 算术平均 ipTM 为 0.8660，尚未达到 0.90。
- 两条Accepted_AF3均为25/25模型两类界面各实际接触≥3个突变位点；OVA-C4-073/AI-MULTIHIT-11的突变—预测界面平均复现比例分别为0.7647/0.7060，ProLIF非纯VdW设计作用位点分别为6/7。

## 权威结果文件

- 旧 Final20 单一总表：`experiments/e151_final20_screen_v1/unified/final20_failclosed_screen.tsv`
- 旧 Final20 阶段计数：`experiments/e151_final20_screen_v1/unified/screen_counts.tsv`
- 当前Accepted_AF3排名CIF、全序列FASTA、含设计界面一致性指标的单一总表、摘要及文件数审计：`deliverables/current_accepted_af3_ranked_20260903.tgz`（新25模型命中或v1.1状态变化后自动原子刷新）
- `OVA-C4-073` 25 模型 AF3 重算：`experiments/e165_c4_iptm90_073_multiseed/final25_screen/af3_candidate_recomputed_summary.tsv`
- `OVA-C4-105` 25 模型失败重算：`experiments/e159_c4_iptm90_105_multiseed/final25_screen/af3_candidate_recomputed_summary.tsv`
- E162 完整 seed1 表：`experiments/e162_c4_105_reversion_ai_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv`
- E174 分片 seed1 表：`experiments/e174_c4_073_multihit_ai_af3_seed1/partial_s*/seed1_strict_summary.tsv`
- E174 `AI-MULTIHIT-01` 25模型严格汇总：`experiments/e183_c4_multihit01_multiseed/full25_summary.tsv`
- E174 `AI-MULTIHIT-11` 25模型严格重算：`experiments/e186_c4_multihit11_multiseed/full25_screen/af3_candidate_recomputed_summary.tsv`
- E174 `AI-MULTIHIT-11` AF3-medoid ProLIF：`experiments/e186_c4_multihit11_multiseed/prolif_af3_medoid/`
- E174 `AI-MULTIHIT-11` Protenix-v2五seed失败汇总：`experiments/e190_c4_multihit11_screen_v1/protenix_summary/protenix_candidate_summary.tsv`
- E174 `AI-MULTIHIT-11` Rosetta 200次失败汇总：`experiments/e190_c4_multihit11_screen_v1/rosetta/ova_body_c4_multihit11_ms_AF3_rosetta_summary.tsv`
- E174 `AI-MULTIHIT-11` OpenDDE十seed失败汇总：`experiments/e203_c4_multihit11_opendde_msa1024_v2/summary/opendde_candidate_summary.tsv`
- E174 `AI-MULTIHIT-11` 当前fail-closed单表：`experiments/e190_c4_multihit11_screen_v1/unified_current_after_opendde/final20_failclosed_screen.tsv`
- `OVA-C4-073` OpenDDE十seed失败汇总：`experiments/e214_c4_073_opendde_v1/summary/opendde_candidate_summary.tsv`
- `OVA-C4-073` Protenix-v2五seed失败汇总：`experiments/e213_c4_073_screen_v1/protenix_summary/protenix_candidate_summary.tsv`
- `OVA-C4-073` 最终fail-closed单表：`experiments/e213_c4_073_screen_v1/unified_final/final20_failclosed_screen.tsv`
- `OVA-C4-073` stoichiometry通过汇总：`experiments/e213_c4_073_screen_v1/stoichiometry/stoichiometry_candidate_summary.tsv`
- `OVA-C4-073` Rosetta 200次失败汇总：`experiments/e213_c4_073_screen_v1/rosetta/ova_body_c4_073_ms_AF3_rosetta_summary.tsv`
- `OVA-C4-073` v1.0报告：`experiments/e213_c4_073_screen_v1/report/筛选总结.md`
- E174 `AI-MULTIHIT-03` 25模型失败汇总：`experiments/e187_c4_multihit03_multiseed/full25_summary.tsv`
- E174 `AI-MULTIHIT-06` 独立seed2失败汇总：`experiments/e188_c4_multihit06_multiseed/seed2_strict/seed1_strict_summary.tsv`
- E174 `AI-MULTIHIT-15` 25模型失败汇总：`experiments/e189_c4_multihit15_multiseed/full25_summary.tsv`
- E177 32条严格seed1汇总：`experiments/e177_c4_tied_proteinmpnn_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv`
- E177候选02/23独立seed2失败：`experiments/e195_c4_e177_02_multiseed/seed2_strict/seed1_strict_summary.tsv`、`experiments/e196_c4_e177_23_multiseed/seed2_strict/seed1_strict_summary.tsv`
- E178完整16条严格seed1汇总：`experiments/e178_c4_073_af3ensemble_ai_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv`
- E178候选04/15独立多seed复核：`experiments/e200_c4_e178_04_multiseed/`、`experiments/e201_c4_e178_15_multiseed/`
- E178候选15局部clash定向tied-AF2修复：`experiments/e204_c4_276_strictensemble_tied_design/`
- E204/E205联合冻结E212 Top32：`experiments/e204_c4_276_strictensemble_tied_design/candidate_pool_joint_e205_novel/`、`experiments/e212_c4_276fix_joint_ai_af3_seed1/`
- OVA-C4-073平衡/高增益tied-AF2连续设计：`experiments/e215_c4_073_balanced_tied_design/`、`experiments/e217_c4_073_highgain_tied_design/`
- OVA-C4-073一recycle设计与E220 AF3链：`experiments/e219_c4_073_recycle1_tied_design/`、`experiments/e220_c4_073_recycle1_ai_af3_seed1/`
- AI11 AF3–Protenix跨模型设计与联合MPNN：`experiments/e221_c4_ai11_af3_protenix_crossmodel_tied_design/`、`experiments/e223_c4_ai11_crossmodel_tied_mpnn/`
- OVA073 Rosetta弱界面极性偏置ProteinMPNN：`experiments/e225_c4_073_rosetta_polar_tied_mpnn/`
- 跨模型主动学习v2、候选权衡图与中卫tied-AF2：`experiments/e234_crossmodel_active_learning_v2/`、`experiments/e235_c4_crossmodel_active_tied_design/`
- 跨OVA073/AI11骨架、弱界面自适应梯度tied-AF2及E244 AF3链：`experiments/e243_c4_crossensemble_active_tied_design/`、`experiments/e244_c4_crossensemble_active_ai_af3_seed1/`
- E223/E225快速MPNN库及E246/E247完整六条件重采样、E248联合AF3链：`experiments/e223_c4_ai11_crossmodel_tied_mpnn/`、`experiments/e225_c4_073_rosetta_polar_tied_mpnn/`、`experiments/e246_c4_ai11_crossmodel_tied_mpnn_full/`、`experiments/e247_c4_073_rosetta_polar_tied_mpnn_full/`、`experiments/e248_c4_fullcondition_mpnn_joint_ai_af3_seed1/`
- 二阶突变主动学习、E251跨骨架设计与优先E252 AF3：`experiments/e250_c4_epistasis_active_learning/`、`experiments/e251_c4_epistasis_seeded_tied_design/`、`experiments/e252_c4_epistasis_ai_af3_seed1/`
- 旧slim archive轨迹恢复、414条增量库与E257 AF3：`experiments/e256_c4_recovered_trajectory_library/`、`experiments/e257_c4_recovered_trajectory_ai_af3_seed1/`
- E235完整轨迹、E236代理优选AF3及E259直接双状态AF2救援AF3：`experiments/e235_c4_crossmodel_active_tied_design/`、`experiments/e236_c4_crossmodel_active_ai_af3_seed1/`、`experiments/e259_c4_direct_tied_af2_rescue_af3_seed1/`
- E243最佳离散状态跨四套AF3骨架续跑及E262自动AF3链：`experiments/e261_c4_e243_crossbackbone_continuation_tied_design/`、`experiments/e262_c4_e243_crossbackbone_continuation_af3_seed1/`
- AI11缺失计量补测与跨骨架共识压缩E265/E266：`experiments/e264_c4_multihit11_stoichiometry_v1/`、`experiments/e265_c4_crossbackbone_consensus_tied_design/`、`experiments/e266_c4_crossbackbone_consensus_af3_seed1/`
- E251近离散软盆地稳定化及E269无模板AF3链：`experiments/e268_c4_e251_softbasin_stabilization_tied_design/`、`experiments/e269_c4_e251_softbasin_stabilization_af3_seed1/`
- 八骨架精确tied-ProteinMPNN重评分及E272无模板AF3链：`experiments/e271_c4_multibackbone_exact_mpnn_score/`、`experiments/e272_c4_multibackbone_exact_mpnn_af3_seed1/`
- E261近离散盆地独立稳定化及E275/E276无模板AF3链：`experiments/e274_c4_e261_softbasin_stabilization_tied_design/`、`experiments/e275_c4_e261_softbasin_stabilization_af3_seed1/`、`experiments/e276_e275_top4_af3_multiseed/`
- E168前四shard严格seed1表：`experiments/e168_c4_073_reversion_ai_af3_seed1/seed1_strict_s*/`
- E168五条独立多seed复核：`experiments/e206_c4_e168_rev002_multiseed/`、`e207`、`e208`、`e209`、`e211`
- 中卫v7新增序列及AI重排：`experiments/e169_c4_073_af3ensemble_tied_design/candidate_pool_v7_novel/`
- 中卫v8新增序列及AI重排：`experiments/e169_c4_073_af3ensemble_tied_design/candidate_pool_v8_novel/`
- 中卫AI11 E193-v1新增序列及AI重排：`experiments/e193_c4_multihit11_af3ensemble_tied_design/candidate_pool_ai11v1_novel/`
- 中卫AI11 E193-v2新增序列：`experiments/e193_c4_multihit11_af3ensemble_tied_design/candidate_pool_ai11v2_novel/`
- 中卫AI11 E193-v3新增序列及加权AF3代理Top32：`experiments/e193_c4_multihit11_af3ensemble_tied_design/candidate_pool_ai11v3_novel/`
- 中卫v6新增序列及AI重排：`experiments/e169_c4_073_af3ensemble_tied_design/candidate_pool_v6_novel/`
- 逐实验假设、失败原因和保留决定：`research_log.md`

## 计算审计备注

### 已完成完整后筛样本的逐阶段汇总

- 统一逐样本结果：`COMPLETED_POSTSCREEN_RESULTS.tsv`
- 可读版报告：`COMPLETED_POSTSCREEN_RESULTS.md`
- 该报告包含 E151 Final20 的20条历史记录，以及 `AI-MULTIHIT-11`、`OVA-C4-073`、`ova_body_c4_06_ms25_r04` 三条单独完成完整后筛的候选；全部 required 通过数仍为0。

## 2026-09-03 · 当前实时快照（覆盖前述历史批次计数）

- **Accepted_AF3：23 条**。最新增加来自 E316 的 `c02/c07/c28/c21` 和 E317 的 `c30/c09/c16/c06/c29`；每条均有无模板 AF3 五个独立 seed、每个 5 samples，共 25 个模型的完整记录，并通过 25/25 零跨链重原子 clash、完整 C4、天然链内二硫键、零亚基间二硫键和双设计界面门槛。
- 当前完整 25 模型均值最高为 `ova_body_c4_30_ms25_r02`：16 个突变、表面突变比例 0.875、平均/最低 ipTM `0.8624/0.8400`。因此“接近或达到 0.90”仍是未完成的优化目标，不能把当前最高值写成 0.90。
- E317 的 3 条未晋级候选（`c14/c15/c32`）保留为几何不稳定负标签；它们不进入 Accepted_AF3 或 v1.0 正向候选池。
- **Accepted_v1.1：目前仍为 0**。23 条 Accepted_AF3 已进入统一 v1.1 后筛队列；ProLIF、Protenix-v2、化学计量和 Rosetta 的 required 项继续 fail-closed 统计，OpenDDE仅作可选诊断。
- 权威实时注册表：`experiments/current_accepted_af3_registry/accepted_af3_registry.tsv`；E317完整汇总：`experiments/e317_e314_rank13_20_af3_multiseed/full25/full25_summary.tsv`；当前排名交付包由 `af3_pipeline/refresh_current_accepted_delivery.sh` 自动刷新。

- 火山 `ky-20260901-558`–`565` 因错误读取旧母本已全部停止，明确排除于所有统计。
- 修复后的 `ky-20260901-566`–`573` 使用 `OVA-C4-073-BB1`–`BB4`，仅在生成完成并通过统一过滤后计入候选池。
## 2026-09-02 E277 / E278：E268完全离散高分状态跨骨架续跑

- E265、E268、E271、E274已完成并回收；其中E268出现当前最强的完全离散生成状态：19个突变、表面比例1.00、tied-AF2两状态ipTM 0.762/0.745、最低pLDDT 0.911。
- E274未保留E261连续态优势，当前不继续该盆地；E271多骨架ProteinMPNN结果只用于E272排序，不计作AF3证据。
- E277任务名为`ky-20260902-1320`至`ky-20260902-1327`，八任务均0退出并完整回收；20突变硬预算、四个独立AF3骨架上下文、每个两种界面权重。最佳单任务hard状态两界面tied-AF2约0.742/0.739，但跨骨架一致性仍差。
- 去除历史已测序列后获得110条有效唯一新序列；AF3-Hamming代理预测仅约0.201–0.349。按代理、距离和多样性冻结21条E278候选，等待火山无模板AF3 seed1；Top4再补齐5 seeds×5 samples。AF3任务号预留1328–1351。
- 当前验收不变：Accepted_AF3=2，Accepted_v1.0=0；E277的tied-AF2和代理值均不计作AF3证据。
## 2026-09-02 E280 / E281：四骨架八状态联合tied-AF2

- 动机：E277初步显示E268单骨架hard高分未稳定迁移；多条轨迹的一侧界面降至0.14–0.27。
- 方法：同一条四链同源tied序列在每次更新中同时通过4套独立AF3骨架的AC/AD共8个二聚状态，所有梯度联合，当前最弱状态梯度×2.5；折叠、突变预算和表面比例同时约束。
- 输入：8条唯一锚序列，来自OVA-C4-073、AI-MULTIHIT-11、E268/E265完全离散状态及E271八骨架ProteinMPNN Top候选；14–20个突变、表面比例0.80–1.00、天然Cys/SIINFEKL/糖基化保持。
- 链映射审计发现第4套骨架的`AD`并非同一物理界面：相对参考第二界面的5.3 Å接触图Jaccard为0，伴侣链相对取向偏差约61.1 Å；同模型`AB`才匹配，Jaccard 0.821。
- `ky-20260902-1352`–`1359`已在首次迭代后全部停止，E280无候选、永久排除于统计；E281取消。该结果是输入链标签错误，不是序列AF3失败。

## 2026-09-02 E283 / E284：接触图重映射后的四骨架八状态联合设计

- 对每套AF3四聚体全部6个链对与参考两类界面做正反向5.3 Å接触图Jaccard匹配，冻结映射`AC/AD、BD/AD、AC/AD、AC/AB`；8个状态的最小Jaccard为0.821，全部超过0.70硬门槛。
- 重新生成8个二聚目标、restraint与36位联合设计mask；8条锚序列及20突变硬上限不变。
- 中卫任务`ky-20260902-1384`–`1391`（Slurm 153771–153778）已全部运行。完成后冻结E284；火山严格seed1预留1392–1399、Top4五seed复核预留1400–1415。
- 火山顺序调整为E252→E284→E287→E257→E259→E236→E254及其余旧批；低历史代理E278仍在E275之后。只有无模板AF3通过后才进入Accepted_AF3与v1.0。

## 2026-09-02 E286 / E287：底部2–3状态稳健续跑

- E283完整归档后自动从hard=1、≤20突变、表面比例≥0.80的精确轨迹中，按8状态最小/平均tied-AF2 ipTM及最低pLDDT选择8个多样锚点。
- 前4任务使用纯表面mask并同时增强当前最低2状态梯度2.2倍；后4任务同时增强最低3状态1.8倍，并开放在四套第一界面中全部复现的热点334/336/338。热点突变仍受核心惩罚和表面比例≥0.80硬门槛。

## Rosetta v1.1 主设计界面口径（2026-09-03）

针对四聚体对角/背景链对可能没有氢键的问题，Rosetta 的 SC、氢键、
dG/dSASA 和埋藏未满足极性原子密度门槛改为只在 AF3 代表结构中设计残基覆盖数最高的一类
主界面的两个对称副本上计算；全部有效链对仍保留用于四聚体 clash、穿插、连通拓扑和逐对诊断。
三条已有200次FastRelax集合按此口径重算后仍未通过：AI-MULTIHIT-11 主界面最弱
SC/H-bond/unsat = 0.687/0/6.40，OVA-C4-073 = 0.673/0/3.71，
c06 = 0.682/0/1.73。也就是说，排除其他界面类型没有改变当前三条候选的Rosetta总通过结论；
三条主设计界面的最弱氢键仍为0。c06的主界面unsat已通过≤2门槛，目前Rosetta只剩氢键门槛失败。
机器可读结果见各候选目录下的 `rosetta_primary_v11/*_rosetta_summary.tsv`，实现脚本为
`af3_pipeline/recompute_rosetta_primary_interface_metrics.py`。
- 中卫任务预留1416–1423，E287火山严格seed1预留1424–1431、Top4五seed复核预留1432–1447；自动启动观察器已运行。当前只是预注册，无候选或AF3结果。

## 2026-09-02 E252火山启动状态

- c20250508空闲GPU达到12后，E252的32条经验AI候选按8个shard启动无模板AF3 seed1。
- 首次任务`ky-20260903-1056`–`1063`因共享磁盘配额在写首个data JSON时全部失败，没有产生AF3模型，不计入候选评价。
- 清理9.084 GiB流水线规则本应删除的旧完整置信度JSON后，以`ky-20260903-1448`–`1455`重提；32/32候选、160/160模型完成。严格seed1命中5条：09/27/30/26/23，平均ipTM依次0.854/0.850/0.850/0.818/0.802，均5/5零clash、完整C4、双设计界面、零链间二硫键。
- Top4候选09/30/27/26已冻结E253；seed2任务`ky-20260903-1064/1068/1072/1076`运行中。当前只是seed1命中，Accepted计数不变。

## 2026-09-02 E289：AF3碰撞定向修复

- E252早期完成的8条均保持5/5完整C4和双设计界面，seed1平均ipTM 0.804–0.850，但0/8达到5/5零clash；三条仅差一个sample。
- 新脚本逐原子定位碰撞，当前复现于95/96–348、187–345、340–345、354–354等设计相关链间对。E289将在E252全部32条结束后，仅对实际碰撞中的已突变位点使用P3-13R回退或Accepted_AF3等位基因做1–3点修复，并将新AF3标签加入ipTM/零clash联合机器学习排序。
- E289从17个有碰撞父本生成127条唯一修复，加入E252新标签后的ipTM/零clash模型交叉验证Spearman为0.770/0.771；排除1548条历史已测序列后冻结Top32。正式序列14–20突变、表面比例0.80–1.00，天然Cys和SIINFEKL全部完整。预留seed1任务1456–1463、Top4多seed任务1464–1479；火山唯一链为E252/E253→E289→E284→E287，当前尚无E289 AF3结果。

## 2026-09-03 E311 / E312：独立四骨架续设计接入无模板AF3

- E311 中卫 8×A800 八状态 tied-AF2 续设计已 8/8 完成并回传；40 条有效轨迹候选均为 13–20 个突变、表面突变比例 0.8667–1.0。
- 结合历史 AF3 标签的加权 Hamming 排序并要求最小距离3后，E312 冻结 11 条独立候选（13–19 个突变、表面比例 0.875–1.0），配置为无模板 AF3 seed1×5 samples、8 分片。
- E312 不再等待旧 E295/E287 链路；火山仍只使用 `q-20251224105744-qhrgl`/`c20250508`，并严格要求 free GPU>10。当前 50/56 张占用、free=6，故尚未提交；watcher 持续轮询。
- 这些候选尚未计入 Accepted_AF3，必须完成 seed1 严格门槛后才可扩展 seed2–5，并继续进入 v1.0 后筛。

## 2026-09-03 · c06 OpenDDE后筛完成

- `ova_body_c4_06_ms25_r04` 的中卫 OpenDDE 任务 `ky-20260903-1616`（Slurm 153994）已正常完成，10/10 seed 模型全部回传。
- OpenDDE 平均/最高 ipTM 为 `0.4124/0.4442`，几何可选代表 `0/10`，因此 OpenDDE required gate 失败；该候选仍不是 `Accepted_v1.0`。
- 回传汇总曾因旧 AF3 表缺少 `topology_cluster`、`min_incident_iptm` 字段中止，已修复兼容性并重跑；修复不改变任何 ipTM、clash、DockQ 或 required 门槛。
- 当前注册表仍为 `Accepted_AF3=23`、`Accepted_v1.1=0`，待筛条目由22条降为21条；OpenDDE不再作为required；交付包已原子刷新。

## 2026-09-03 · E318/E321部分接力状态

- E318的7个中卫shard已完成并回收，部分候选池有37条有效唯一序列（13–20个突变、表面突变比例0.867–1.0）；第8个shard仍在运行，因此这些序列尚未进入AF3注册表。
- 为避免释放的A800闲置，已从部分池按历史AF3标签和Hamming距离冻结8个E321锚点，启动`ky-20260903-2500`–`2506`（Slurm 154107–154113）；`2507`（154114）等待E318最后任务释放资源。
- E321此前因manifest路径缺少`experiments/`前缀而失败的首轮任务已保留为失败审计；路径已修复并重新提交。E319仍严格等待完整E318八shard归档，E321部分池不会替代完整AF3证据。
