#!/usr/bin/env bash
# Validate the AI-ranked clean-chassis allele bridge after E301 multiseed.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
PREVIOUS="$PROJECT/experiments/e302_e301_top4_af3_multiseed"
EXPERIMENT="$PROJECT/experiments/e303_c4_cleanbridge_allele_ai_af3_seed1"
PY=/root/400083/alphafold3/repo/.venv/bin/python
cd "$PROJECT"
while [[ ! -f "$PREVIOUS/MULTISEED_VALIDATION_COMPLETE" || \
         ! -f "$EXPERIMENT/CANDIDATES_FROZEN" ]]; do sleep 30; done
expected_models=$(( $(tail -n +2 "$EXPERIMENT/manifest.tsv" | wc -l) * 5 ))
while [[ ! -f "$EXPERIMENT/VOLC_AF3_SHARDS_SUBMITTED" ]]; do
  if [[ -f "$EXPERIMENT/VOLC_AF3_PARTIAL7_SUBMITTED" ]]; then
    "$PY" af3_pipeline/submit_af3_volc_if_capacity.py \
      --experiment "$EXPERIMENT" --job-number-start 1760 \
      --shards 7 --expected 8 --record-file volc_af3_jobs.tsv \
      --marker VOLC_AF3_SHARDS_SUBMITTED --priority
  else
    "$PY" af3_pipeline/submit_af3_volc_if_capacity.py \
      --experiment "$EXPERIMENT" --job-number-start 1760 --priority
  fi
  [[ -f "$EXPERIMENT/VOLC_AF3_SHARDS_SUBMITTED" ]] && break
  sleep 30
done
while true; do
  complete=$(find "$EXPERIMENT" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  [[ "$complete" -eq "$expected_models" ]] && break
  sleep 30
done
"$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
  --experiment "$EXPERIMENT" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
  --out "$EXPERIMENT/seed1_strict_final"
touch "$EXPERIMENT/SEED1_STRICT_RESULTS_COMPLETE"
bash af3_pipeline/run_seed1_hits_multiseed_volc.sh \
  "$EXPERIMENT" "$PROJECT/experiments/e304_e303_top4_af3_multiseed" 1768 4
touch "$EXPERIMENT/SEED1_STRICT_SCREEN_COMPLETE"
echo "E303 clean-bridge AF3 screen complete"
