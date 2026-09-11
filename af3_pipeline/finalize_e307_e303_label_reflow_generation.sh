#!/usr/bin/env bash
# Re-rank the untested clean-chassis allele space after adding all E303 AF3 labels.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
REGISTRY="$PROJECT/experiments/current_accepted_af3_registry/accepted_af3_registry.tsv"
WORK="$PROJECT/experiments/e307_c4_e303_label_reflow_generation"
OUT="$PROJECT/experiments/e307_c4_e303_label_reflow_af3_seed1"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"
if [[ -f "$OUT/CANDIDATES_FROZEN" ]]; then echo "E307 already frozen"; exit 0; fi
mkdir -p "$WORK"

"$PY" af3_pipeline/make_multihit_allele_library.py \
  --manifest "$REGISTRY" \
  --system ova_body_c4_05_ms25_r01 \
  --system ova_body_c4_01_ms25_r02 \
  --system ova_body_c4_073_ms \
  --system ova_body_c4_06_ms25_r04 \
  --system ova_body_c4_multihit11_ms \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --surface-positions experiments/e283_c4_chainmapped_eightstate_tied_design/surface_positions.txt \
  --exclude-glob 'experiments/e*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/*/manifest.tsv' \
  --out "$WORK" --max-library 4096 --name-prefix E307-REFLOW

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
  "${ranking_args[@]}" --library "$WORK/allele_library.tsv" \
  --out "$WORK/hamming_ai_ranking.tsv" --shortlist 32 --minimum-distance 2 \
  --allow-short-shortlist --duplicate-label-policy max_models_mean --weighted-iptm
"$PY" af3_pipeline/prepare_ova_body_af3.py \
  --fasta "$WORK/hamming_ai_ranking_shortlist.fasta" \
  --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
  --out "$OUT" --top 32 --seeds 1 --seed-start 1 --seed-offset 0 \
  --shards 8 --n-chains 4
touch "$OUT/CANDIDATES_FROZEN"
echo "E307 E303-label reflow candidates frozen"
