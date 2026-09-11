#!/usr/bin/env bash
# Freeze E312 AF3 candidates from the continued E311 Zhongwei trajectories.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E311="$PROJECT/experiments/e311_c4_e294_detbest_continuation_tied_design"
E312="$PROJECT/experiments/e312_c4_e311_detbest_af3_seed1"
POOL="$E311/candidate_pool_novel"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"
if [[ -f "$E312/CANDIDATES_FROZEN" ]]; then echo "E312 already frozen"; exit 0; fi
while [[ ! -f "$E311/REMOTE_ARCHIVE_RECOVERED" ]]; do sleep 30; done

"$PY" af3_pipeline/collect_iptm90_crossmodel_candidates.py \
  --experiment "$E311" --run-glob 'run_bottomk_*' \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
  --exclude-glob 'experiments/e*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/*/manifest.tsv' \
  --exclude-glob 'deliverables/**/*.fasta' \
  --exclude-ignore-path "$E312/manifest.tsv" \
  --out "$POOL" --top-af2 256 --top-surrogate 0 --top-mpnn 0 \
  --minimum-distance 2 --min-mutations 1 --max-mutations 20

mapfile -t rankings < <(awk '
  /^rankings=\(/{inside=1; next}
  inside && /^\)/{exit}
  inside {gsub(/^[[:space:]]+|[[:space:]]+$/, ""); if (length) print}
' af3_pipeline/finalize_e217_generation.sh)
ranking_args=(); for file in "${rankings[@]}"; do ranking_args+=(--ranking "$file"); done
while IFS= read -r file; do ranking_args+=(--ranking "$file"); done < <(
  find experiments -type f \( \
    -path '*/seed1_strict_final/seed1_strict_summary.tsv' -o \
    -path '*/full25/full25_summary.tsv' \
  \) | sort
)
"$PY" af3_pipeline/score_c4_library_hamming_ai.py \
  "${ranking_args[@]}" --library "$POOL/generation_shortlist.tsv" \
  --out "$POOL/hamming_ai_ranking_current_weighted.tsv" \
  --shortlist 32 --minimum-distance 3 --allow-short-shortlist \
  --duplicate-label-policy max_models_mean --weighted-iptm \
  --generation-rank-weight 0.10
"$PY" af3_pipeline/prepare_ova_body_af3.py \
  --fasta "$POOL/hamming_ai_ranking_current_weighted_shortlist.fasta" \
  --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
  --out "$E312" --top 32 --seeds 1 --seed-start 1 --seed-offset 0 \
  --shards 8 --n-chains 4
touch "$E312/CANDIDATES_FROZEN"
echo "E311 trajectories frozen as E312"
