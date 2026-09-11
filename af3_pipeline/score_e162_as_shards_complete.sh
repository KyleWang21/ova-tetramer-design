#!/usr/bin/env bash
# Materialize frozen per-shard AF3 geometry gates as E162 finishes.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
EXP="$PROJECT/experiments/e162_c4_105_reversion_ai_af3_seed1"
PY=/root/400083/alphafold3/repo/.venv/bin/python
cd "$PROJECT"

for shard in 0 1 2 3 4 5 6 7; do
  while true; do
    found=$(find "$EXP/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
    [[ "$found" -eq 10 ]] && break
    sleep 30
  done
  "$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
    --experiment "$EXP" \
    --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
    --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
    --shard "$shard" --out "$EXP/partial_s$shard"
done

"$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
  --experiment "$EXP" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
  --out "$EXP/seed1_strict_final"
