#!/usr/bin/env bash
# Filter and AF3-label-rank the OVA-C4-073 highgain-v10 tied-AF2 library.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E217="$PROJECT/experiments/e217_c4_073_highgain_tied_design"
E218="$PROJECT/experiments/e218_c4_073_highgain_v10_ai_af3_seed1"
POOL="$E217/candidate_pool_novel"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"
if [[ -f "$E218/CANDIDATES_FROZEN" ]]; then
  echo "E218 already frozen; nothing to do"
  exit 0
fi

while [[ ! -f "$E217/REMOTE_ARCHIVE_RECOVERED" ]]; do
  sleep 30
done
"$PY" af3_pipeline/collect_iptm90_crossmodel_candidates.py \
  --experiment "$E217" --run-glob 'run_073v10_*' \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
  --exclude-glob 'experiments/e*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/*/manifest.tsv' \
  --exclude-glob 'experiments/e193_c4_multihit11_af3ensemble_tied_design/candidate_pool_ai11v*_novel/*.fasta' \
  --exclude-glob 'experiments/e205_c4_276_strictensemble_tied_mpnn/candidate_pool_novel/*.fasta' \
  --exclude-glob 'experiments/e215_c4_073_balanced_tied_design/candidate_pool_novel/*.fasta' \
  --exclude-glob 'deliverables/**/*.fasta' \
  --exclude-ignore-path "$E218/manifest.tsv" \
  --out "$POOL" --top-af2 256 --top-surrogate 0 --top-mpnn 0 \
  --minimum-distance 2 --min-mutations 1 --max-mutations 24

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
for file in "${rankings[@]}"; do
  ranking_args+=(--ranking "$file")
done
RANKED="$POOL/hamming_ai_ranking_current_weighted.tsv"
"$PY" af3_pipeline/score_c4_library_hamming_ai.py \
  "${ranking_args[@]}" --library "$POOL/generation_shortlist.tsv" \
  --out "$RANKED" --shortlist 32 --minimum-distance 3 \
  --allow-short-shortlist --duplicate-label-policy max_models_mean --weighted-iptm

"$PY" af3_pipeline/prepare_ova_body_af3.py \
  --fasta "$POOL/hamming_ai_ranking_current_weighted_shortlist.fasta" \
  --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
  --out "$E218" --top 32 --seeds 1 --seed-start 1 --seed-offset 0 \
  --shards 8 --n-chains 4
cp af3_pipeline/e218_jobs.tsv.template "$E218/jobs.tsv"
touch "$E218/CANDIDATES_FROZEN"
echo "E217 filtered, ranked, and frozen as E218 Top32"
