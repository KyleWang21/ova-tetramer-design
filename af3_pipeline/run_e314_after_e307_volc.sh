#!/usr/bin/env bash
# Use the completed E307 AF3 labels for one more finite-library Bayesian-style
# acquisition round before returning Fire GPUs to older exploratory branches.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
PREVIOUS="$PROJECT/experiments/e307_c4_e303_label_reflow_af3_seed1"
WORK="$PROJECT/experiments/e314_c4_e307_epistasis_reflow_generation"
EXPERIMENT="$PROJECT/experiments/e314_c4_e307_epistasis_reflow_af3_seed1"
REGISTRY="$PROJECT/experiments/current_accepted_af3_registry/accepted_af3_registry.tsv"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
AF3_PY=/root/400083/alphafold3/repo/.venv/bin/python
cd "$PROJECT"

while [[ ! -f "$PREVIOUS/SEED1_STRICT_SCREEN_COMPLETE" ]]; do sleep 30; done

if [[ ! -f "$EXPERIMENT/CANDIDATES_FROZEN" ]]; then
  mkdir -p "$WORK"
  "$PY" af3_pipeline/make_multihit_allele_library.py \
    --manifest "$REGISTRY" \
    --system ova_body_c4_22_ms25_r01 \
    --system ova_body_c4_12_ms25_r03 \
    --system ova_body_c4_10_ms25_r03 \
    --system ova_body_c4_073_ms \
    --system ova_body_c4_06_ms25_r04 \
    --system ova_body_c4_multihit11_ms \
    --system ova_body_c4_13_ms25_r06 \
    --system ova_body_c4_05_ms25_r01 \
    --system ova_body_c4_01_ms25_r02 \
    --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
    --surface-positions experiments/e283_c4_chainmapped_eightstate_tied_design/surface_positions.txt \
    --exclude-glob 'experiments/e*/manifest.tsv' \
    --exclude-glob 'experiments/e*/*/manifest.tsv' \
    --exclude-glob 'experiments/e*/*/*/manifest.tsv' \
    --include-reference-position 354 \
    --out "$WORK" --max-library 4096 --name-prefix E314-EPIREFLOW

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
  "$PY" af3_pipeline/rerank_surrogate_epistasis.py \
    "${ranking_args[@]}" --library "$WORK/hamming_ai_ranking.tsv" \
    --out "$WORK" --top 32 --minimum-distance 2 --library-manifold-only
  "$PY" af3_pipeline/prepare_ova_body_af3.py \
    --fasta "$WORK/epistasis_reranked_candidates.fasta" \
    --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
    --out "$EXPERIMENT" --top 32 --seeds 1 --seed-start 1 --seed-offset 0 \
    --shards 8 --n-chains 4
  touch "$EXPERIMENT/CANDIDATES_FROZEN"
fi

expected_models=$(( $(tail -n +2 "$EXPERIMENT/manifest.tsv" | wc -l) * 5 ))
while [[ ! -f "$EXPERIMENT/VOLC_AF3_SHARDS_SUBMITTED" ]]; do
  "$AF3_PY" af3_pipeline/submit_af3_volc_if_capacity.py \
    --experiment "$EXPERIMENT" --job-number-start 2240
  [[ -f "$EXPERIMENT/VOLC_AF3_SHARDS_SUBMITTED" ]] && break
  sleep 30
done
while true; do
  complete=$(find "$EXPERIMENT" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  [[ "$complete" -eq "$expected_models" ]] && break
  sleep 30
done
"$AF3_PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
  --experiment "$EXPERIMENT" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
  --out "$EXPERIMENT/seed1_strict_final"
touch "$EXPERIMENT/SEED1_STRICT_RESULTS_COMPLETE"
bash af3_pipeline/run_seed1_hits_multiseed_volc.sh \
  "$EXPERIMENT" "$PROJECT/experiments/e315_e314_top4_af3_multiseed" 2248 4
touch "$EXPERIMENT/SEED1_STRICT_SCREEN_COMPLETE"
echo "E314 E307-label epistasis reflow AF3 screen complete"
