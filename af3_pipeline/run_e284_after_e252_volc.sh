#!/usr/bin/env bash
# Prioritize corrected chain-mapped E284 after all E289 seed1 hits finish multiseed validation.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
PREVIOUS="$PROJECT/experiments/e291_e289_remaining5_af3_multiseed"
EXPERIMENT="$PROJECT/experiments/e284_c4_chainmapped_eightstate_af3_seed1"
PY=/root/400083/alphafold3/repo/.venv/bin/python
cd "$PROJECT"
while [[ ! -f "$PREVIOUS/MULTISEED_VALIDATION_COMPLETE" || \
         ! -f "$EXPERIMENT/CANDIDATES_FROZEN" ]]; do sleep 30; done
expected_models=$(( $(tail -n +2 "$EXPERIMENT/manifest.tsv" | wc -l) * 5 ))
while [[ ! -f "$EXPERIMENT/VOLC_AF3_SHARDS_SUBMITTED" ]]; do
  "$PY" af3_pipeline/submit_af3_volc_if_capacity.py \
    --experiment "$EXPERIMENT" --job-number-start 1392
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
  "$EXPERIMENT" "$PROJECT/experiments/e285_e284_top4_af3_multiseed" 1400 4
touch "$EXPERIMENT/SEED1_STRICT_SCREEN_COMPLETE"
echo "E284 corrected chain-mapped eight-state AF3 screen complete after E252"
