#!/usr/bin/env bash
# Recover high-quality near-discrete trajectory states omitted by old slim archives.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
E257="$PROJECT/experiments/e257_c4_recovered_trajectory_ai_af3_seed1"
MERGED="$PROJECT/experiments/e256_c4_recovered_trajectory_library"
cd "$PROJECT"
if [[ -f "$E257/CANDIDATES_FROZEN" ]]; then
  echo "E257 already frozen"
  exit 0
fi

common=(
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif
  --exclude-glob 'experiments/e*/manifest.tsv'
  --exclude-glob 'experiments/e*/*/manifest.tsv'
  --exclude-glob 'experiments/e*/*/*/manifest.tsv'
  --exclude-glob 'deliverables/**/*.fasta'
  --exclude-ignore-path "$E257/manifest.tsv"
  --top-af2 256 --top-surrogate 0 --top-mpnn 0
  --minimum-distance 2 --min-mutations 1 --max-mutations 24
)

declare -a specs=(
  'e217_c4_073_highgain_tied_design|run_073v10_*'
  'e219_c4_073_recycle1_tied_design|run_073r1_*'
  'e221_c4_ai11_af3_protenix_crossmodel_tied_design|run_ai11xm1_*'
)
libraries=()
for spec in "${specs[@]}"; do
  experiment="${spec%%|*}"
  run_glob="${spec#*|}"
  out="$PROJECT/experiments/$experiment/candidate_pool_trajectory_recovered_novel"
  "$PY" af3_pipeline/collect_iptm90_crossmodel_candidates.py \
    --experiment "$PROJECT/experiments/$experiment" --run-glob "$run_glob" \
    --out "$out" "${common[@]}"
  libraries+=(--library "$out/generation_shortlist.tsv")
done

rm -f "$MERGED/merged_recovered_trajectory_library.tsv" \
  "$MERGED/merged_recovered_trajectory_library.fasta"
"$PY" af3_pipeline/merge_recovered_trajectory_libraries.py \
  "${libraries[@]}" --out "$MERGED"

rankings=(
  experiments/e110_lowmut_c4_af3_screen/body_candidate_ranking.tsv
  experiments/e116_c4_iptm80_recomb_af3/body_candidate_ranking.tsv
  experiments/e117_c4_iptm80_af2traj_af3/body_candidate_ranking.tsv
  experiments/e126_c4_iptm80_continue_af3_seed1/body_candidate_ranking.tsv
  experiments/e158_c4_iptm90_joint_surrogate_af3_seed1/seed1_strict_summary.tsv
  experiments/e162_c4_105_reversion_ai_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv
  experiments/e174_c4_073_multihit_ai_af3_seed1/seed1_strict_final_classified/seed1_strict_summary.tsv
  experiments/e177_c4_tied_proteinmpnn_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv
  experiments/e178_c4_073_af3ensemble_ai_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv
  experiments/e168_c4_073_reversion_ai_af3_seed1/seed1_strict_s0/seed1_strict_summary.tsv
  experiments/e168_c4_073_reversion_ai_af3_seed1/seed1_strict_s1/seed1_strict_summary.tsv
  experiments/e168_c4_073_reversion_ai_af3_seed1/seed1_strict_s2/seed1_strict_summary.tsv
  experiments/e168_c4_073_reversion_ai_af3_seed1/seed1_strict_s3/seed1_strict_summary.tsv
  experiments/e195_c4_e177_02_multiseed/seed2_strict/seed1_strict_summary.tsv
  experiments/e196_c4_e177_23_multiseed/seed2_strict/seed1_strict_summary.tsv
  experiments/e200_c4_e178_04_multiseed/seed2_strict/seed1_strict_summary.tsv
  experiments/e201_c4_e178_15_multiseed/seed2_strict/seed1_strict_summary.tsv
  experiments/e159_c4_iptm90_105_multiseed/final25_screen/af3_candidate_recomputed_summary_with_sequence.tsv
  experiments/e165_c4_iptm90_073_multiseed/final25_screen/af3_candidate_recomputed_summary_with_sequence.tsv
  experiments/e183_c4_multihit01_multiseed/full25_summary.tsv
  experiments/e186_c4_multihit11_multiseed/full25_summary.tsv
  experiments/e187_c4_multihit03_multiseed/full25_summary.tsv
  experiments/e189_c4_multihit15_multiseed/full25_summary.tsv
  experiments/e206_c4_e168_rev002_multiseed/full25_strict/combined_25model_summary.tsv
  experiments/e207_c4_e168_rev010_multiseed/seed2_strict/seed1_strict_summary.tsv
  experiments/e208_c4_e168_rev003_multiseed/seed2_strict/seed1_strict_summary.tsv
  experiments/e208_c4_e168_rev003_multiseed/seed3_strict/seed1_strict_summary.tsv
  experiments/e209_c4_e168_rev006_multiseed/seed2_strict/seed1_strict_summary.tsv
  experiments/e211_c4_e168_rev001_multiseed/seed2_strict/seed1_strict_summary.tsv
)
ranking_args=()
for file in "${rankings[@]}"; do ranking_args+=(--ranking "$file"); done
"$PY" af3_pipeline/score_c4_library_hamming_ai.py \
  "${ranking_args[@]}" \
  --library "$MERGED/merged_recovered_trajectory_library.tsv" \
  --out "$MERGED/hamming_ai_ranking.tsv" --shortlist 32 --minimum-distance 3 \
  --allow-short-shortlist --duplicate-label-policy max_models_mean --weighted-iptm

"$PY" af3_pipeline/prepare_ova_body_af3.py \
  --fasta "$MERGED/hamming_ai_ranking_shortlist.fasta" \
  --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
  --out "$E257" --top 32 --seeds 1 --seed-start 1 --seed-offset 0 \
  --shards 8 --n-chains 4
touch "$E257/CANDIDATES_FROZEN"
echo "E257 recovered-trajectory candidates frozen"
