#!/usr/bin/env bash
# Freeze E287 from bottom-k robust chain-mapped trajectories.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E286="$PROJECT/experiments/e286_c4_chainmapped_bottomk_tied_design"
E287="$PROJECT/experiments/e287_c4_chainmapped_bottomk_af3_seed1"
POOL="$E286/candidate_pool_novel"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"
if [[ -f "$E287/CANDIDATES_FROZEN" ]]; then echo "E287 already frozen"; exit 0; fi
while [[ ! -f "$E286/REMOTE_ARCHIVE_RECOVERED" ]]; do sleep 30; done
"$PY" af3_pipeline/collect_iptm90_crossmodel_candidates.py \
  --experiment "$E286" --run-glob 'run_bottomk_*' \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
  --exclude-glob 'experiments/e*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/*/manifest.tsv' \
  --exclude-glob 'experiments/e*/candidate_pool*/*.fasta' \
  --exclude-glob 'experiments/e*/*/candidate_pool*/*.fasta' \
  --exclude-glob 'deliverables/**/*.fasta' \
  --exclude-ignore-path "$E287/manifest.tsv" \
  --out "$POOL" --top-af2 256 --top-surrogate 0 --top-mpnn 0 \
  --minimum-distance 2 --min-mutations 1 --max-mutations 20
mapfile -t rankings < <(awk '
  /^rankings=\(/{inside=1; next}
  inside && /^\)/{exit}
  inside {gsub(/^[[:space:]]+|[[:space:]]+$/, ""); if (length) print}
' af3_pipeline/finalize_e217_generation.sh)
ranking_args=(); for file in "${rankings[@]}"; do ranking_args+=(--ranking "$file"); done
# The legacy ranking list ends at E211.  Add the completed recent AF3
# experiments so E287 is selected from the current measured landscape rather
# than from a stale surrogate.  Missing tables are intentionally skipped: the
# generation watcher may reach this point while a later AF3 batch is running.
recent_rankings=(
  experiments/e252_c4_epistasis_ai_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv
  experiments/e253_e252_top4_af3_multiseed/full25/full25_summary.tsv
  experiments/e284_c4_chainmapped_eightstate_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv
  experiments/e289_c4_e252_clashrepair_ai_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv
  experiments/e290_e289_top4_af3_multiseed/full25/full25_summary.tsv
  experiments/e291_e289_remaining5_af3_multiseed/full25/full25_summary.tsv
  experiments/e292_c4_e283_deterministic_final_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv
  experiments/e297_c4_e284_clashrepair_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv
  experiments/e298_e297_top4_af3_multiseed/full25/full25_summary.tsv
  experiments/e299_c4_e292_clashrepair_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv
  experiments/e300_e299_top4_af3_multiseed/full25/full25_summary.tsv
  experiments/e301_c4_e298_nativeclash_repair_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv
  experiments/e302_e301_top4_af3_multiseed/full25/full25_summary.tsv
  experiments/e303_c4_cleanbridge_allele_ai_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv
  experiments/e304_e303_top4_af3_multiseed/full25/full25_summary.tsv
)
for file in "${recent_rankings[@]}"; do
  [[ -s "$file" ]] && ranking_args+=(--ranking "$file")
done
"$PY" af3_pipeline/score_c4_library_hamming_ai.py \
  "${ranking_args[@]}" --library "$POOL/generation_shortlist.tsv" \
  --out "$POOL/hamming_ai_ranking_current_weighted.tsv" \
  --shortlist 32 --minimum-distance 3 --allow-short-shortlist \
  --duplicate-label-policy max_models_mean --weighted-iptm \
  --generation-rank-weight 0.10
"$PY" af3_pipeline/prepare_ova_body_af3.py \
  --fasta "$POOL/hamming_ai_ranking_current_weighted_shortlist.fasta" \
  --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
  --out "$E287" --top 32 --seeds 1 --seed-start 1 --seed-offset 0 \
  --shards 8 --n-chains 4
touch "$E287/CANDIDATES_FROZEN"
echo "E286 bottom-k robust trajectories frozen as E287"
