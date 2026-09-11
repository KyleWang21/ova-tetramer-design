#!/usr/bin/env bash
# Prioritize the cross-validated epistasis/Hamming ensemble before older AF3 waves.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
EXPERIMENT="$PROJECT/experiments/e252_c4_epistasis_ai_af3_seed1"
PY=/root/400083/alphafold3/repo/.venv/bin/python
cd "$PROJECT"
while [[ ! -f "$EXPERIMENT/CANDIDATES_FROZEN" ]]; do sleep 30; done
expected_models=$(( $(tail -n +2 "$EXPERIMENT/manifest.tsv" | wc -l) * 5 ))
while [[ ! -f "$EXPERIMENT/VOLC_AF3_SHARDS_SUBMITTED" ]]; do
  "$PY" af3_pipeline/submit_af3_volc_if_capacity.py \
    --experiment "$EXPERIMENT" --job-number-start 1448
  [[ -f "$EXPERIMENT/VOLC_AF3_SHARDS_SUBMITTED" ]] && break
  sleep 30
done
while true; do
  if compgen -G "$EXPERIMENT/WORKER_FAILED_s*" >/dev/null; then
    echo "E252 AF3 worker failure marker detected; refusing to wait indefinitely" >&2
    exit 1
  fi
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
  "$EXPERIMENT" "$PROJECT/experiments/e253_e252_top4_af3_multiseed" 1064 4
touch "$EXPERIMENT/SEED1_STRICT_SCREEN_COMPLETE"
echo "E252 epistasis-priority AF3 and promoted seeds2-5 screen complete"
