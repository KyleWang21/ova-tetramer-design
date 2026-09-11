#!/usr/bin/env bash
# Fail-fast multiseed check for the 14-mutation E162 REV036 hit after E182.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E182="$PROJECT/experiments/e182_c4_105_rev88_155_337_356_multiseed"
E180="$PROJECT/experiments/e180_c4_105_rev155_337_356_multiseed"
PY=/root/400083/alphafold3/repo/.venv/bin/python
cd "$PROJECT"

deadline=$((SECONDS + 86400))
while true; do
  [[ -f "$E182/RESEED_DECISION_COMPLETE" ]] && break
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for E182 reseed decision" >&2
    exit 2
  fi
  sleep 30
done

found=$(find "$E180/out_s0" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
if [[ "$found" -ne 5 ]]; then
  bash af3_pipeline/run_worker.sh "$E180" 0
fi
"$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
  --experiment "$E180" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
  --shard 0 --out "$E180/seed2_strict"

passed=$(awk -F '\t' '
  NR==1 {for (i=1; i<=NF; i++) if ($i=="seed1_strict_pass") column=i; next}
  NR==2 {print $column}
' "$E180/seed2_strict/seed1_strict_summary.tsv")
if [[ "$passed" != 1 ]]; then
  echo "E180 failed independent seed2; skipping seeds3-5"
  touch "$E180/RESEED_DECISION_COMPLETE"
  exit 0
fi

for shard in 1 2 3; do
  found=$(find "$E180/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -ne 5 ]]; then
    bash af3_pipeline/run_worker.sh "$E180" "$shard"
  fi
done
touch "$E180/RESEED_DECISION_COMPLETE"
echo "E180 seeds2-5 complete"
