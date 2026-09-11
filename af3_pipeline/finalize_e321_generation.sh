#!/usr/bin/env bash
# E321 archive -> E322 template-free AF3 -> E323 five-seed validation.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E321="$PROJECT/experiments/e321_c4_e318_next_tied_design"
POOL="$E321/candidate_pool_novel"
AF3="$PROJECT/experiments/e322_c4_e321_next_af3_seed1"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
PY_AF3=/root/400083/alphafold3/repo/.venv/bin/python
cd "$PROJECT"

while [[ ! -f "$E321/REMOTE_ARCHIVE_RECOVERED" ]]; do sleep 30; done
if [[ ! -f "$E321/CANDIDATES_FROZEN" ]]; then
  mkdir -p "$POOL"
  "$PY" af3_pipeline/collect_iptm90_crossmodel_candidates.py \
    --experiment "$E321" --run-glob 'run_next_*' \
    --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
    --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
    --exclude-glob 'experiments/e*/manifest.tsv' \
    --exclude-glob 'experiments/e*/*/manifest.tsv' \
    --exclude-glob 'experiments/e*/*/*/manifest.tsv' \
    --exclude-glob 'deliverables/**/*.fasta' \
    --exclude-ignore-path "$E321/manifest.tsv" \
    --out "$POOL" --top-af2 256 --top-surrogate 0 --top-mpnn 0 \
    --minimum-distance 2 --min-mutations 1 --max-mutations 24

  mapfile -t rankings < <(awk '
    /^rankings=\(/{inside=1; next}
    inside && /^\)/{exit}
    inside {gsub(/^[[:space:]]+|[[:space:]]+$/, ""); if (length) print}
  ' af3_pipeline/finalize_e217_generation.sh)
  ranking_args=()
  for file in "${rankings[@]}"; do [[ -f "$file" ]] && ranking_args+=(--ranking "$file"); done
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
  touch "$E321/CANDIDATES_FROZEN"
fi

if [[ ! -f "$AF3/CANDIDATES_FROZEN" ]]; then
  if [[ ! -f "$AF3/manifest.tsv" ]]; then
    "$PY_AF3" af3_pipeline/prepare_ova_body_af3.py \
      --fasta "$POOL/hamming_ai_ranking_current_weighted_shortlist.fasta" \
      --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
      --out "$AF3" --top 32 --seeds 1 --seed-start 1 --seed-offset 0 \
      --shards 8 --n-chains 4
  fi
  [[ -s "$AF3/manifest.tsv" ]] || { echo "E322 manifest missing" >&2; exit 2; }
  touch "$AF3/CANDIDATES_FROZEN"
fi

expected_models=$(( $(tail -n +2 "$AF3/manifest.tsv" | wc -l) * 5 ))
while [[ ! -f "$AF3/VOLC_AF3_SHARDS_SUBMITTED" ]]; do
  "$PY_AF3" af3_pipeline/submit_af3_volc_if_capacity.py \
    --experiment "$AF3" --job-number-start 2600
  [[ -f "$AF3/VOLC_AF3_SHARDS_SUBMITTED" ]] && break
  sleep 30
done
while true; do
  complete=$(find "$AF3" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  [[ "$complete" -eq "$expected_models" ]] && break
  sleep 30
done
"$PY_AF3" af3_pipeline/summarize_iptm90_seed1_strict.py \
  --experiment "$AF3" --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
  --out "$AF3/seed1_strict_final"
touch "$AF3/SEED1_STRICT_RESULTS_COMPLETE"
bash af3_pipeline/run_seed1_hits_multiseed_volc.sh \
  "$AF3" "$PROJECT/experiments/e323_e322_top4_af3_multiseed" 2632 4
touch "$AF3/SEED1_STRICT_SCREEN_COMPLETE"
echo "E321 generation -> E322 AF3 seed1 -> E323 multiseed chain complete"
