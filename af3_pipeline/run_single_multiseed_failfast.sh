#!/usr/bin/env bash
# Run seed2 for one candidate; run seeds3-5 only if seed2 passes strict gates.
set -euo pipefail

EXPERIMENT="${1:?experiment path required}"
PREREQUISITE="${2:?prerequisite marker required}"
PROJECT=/root/400083/ova_p3_tetramer_design
PY=/root/400083/alphafold3/repo/.venv/bin/python
cd "$PROJECT"

deadline=$((SECONDS + 172800))
while [[ ! -f "$PREREQUISITE" ]]; do
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for $PREREQUISITE" >&2
    exit 2
  fi
  sleep 30
done

found=$(find "$EXPERIMENT/out_s0" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
if [[ "$found" -ne 5 ]]; then
  bash af3_pipeline/run_worker.sh "$EXPERIMENT" 0
fi
"$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
  --experiment "$EXPERIMENT" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
  --shard 0 --out "$EXPERIMENT/seed2_strict"

passed=$(awk -F '\t' '
  NR==1 {for (i=1; i<=NF; i++) if ($i=="seed1_strict_pass") column=i; next}
  NR==2 {print $column}
' "$EXPERIMENT/seed2_strict/seed1_strict_summary.tsv")
if [[ "$passed" != 1 ]]; then
  echo "$EXPERIMENT failed independent seed2; skipping seeds3-5"
  touch "$EXPERIMENT/RESEED_DECISION_COMPLETE"
  exit 0
fi

for shard in 1 2 3; do
  found=$(find "$EXPERIMENT/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -ne 5 ]]; then
    bash af3_pipeline/run_worker.sh "$EXPERIMENT" "$shard"
  fi
  seed=$((shard + 2))
  "$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
    --experiment "$EXPERIMENT" \
    --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
    --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
    --shard "$shard" --out "$EXPERIMENT/seed${seed}_strict"
  passed=$(awk -F '\t' '
    NR==1 {for (i=1; i<=NF; i++) if ($i=="seed1_strict_pass") column=i; next}
    NR==2 {print $column}
  ' "$EXPERIMENT/seed${seed}_strict/seed1_strict_summary.tsv")
  if [[ "$passed" != 1 ]]; then
    echo "$EXPERIMENT failed independent seed$seed; skipping later seeds"
    touch "$EXPERIMENT/RESEED_DECISION_COMPLETE"
    exit 0
  fi
done
touch "$EXPERIMENT/RESEED_DECISION_COMPLETE"
echo "$EXPERIMENT seeds2-5 complete"
