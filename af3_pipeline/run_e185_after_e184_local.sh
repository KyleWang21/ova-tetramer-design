#!/usr/bin/env bash
# Fail-fast multiseed check for E174 AI-MULTIHIT-02 after E184.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E184="$PROJECT/experiments/e184_c4_multihit09_multiseed"
E185="$PROJECT/experiments/e185_c4_multihit02_multiseed"
PY=/root/400083/alphafold3/repo/.venv/bin/python
cd "$PROJECT"

deadline=$((SECONDS + 172800))
while true; do
  [[ -f "$E184/RESEED_DECISION_COMPLETE" ]] && break
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for E184 reseed decision" >&2
    exit 2
  fi
  sleep 30
done

found=$(find "$E185/out_s0" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
if [[ "$found" -ne 5 ]]; then
  bash af3_pipeline/run_worker.sh "$E185" 0
fi
"$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
  --experiment "$E185" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
  --shard 0 --out "$E185/seed2_strict"

passed=$(awk -F '\t' '
  NR==1 {for (i=1; i<=NF; i++) if ($i=="seed1_strict_pass") column=i; next}
  NR==2 {print $column}
' "$E185/seed2_strict/seed1_strict_summary.tsv")
if [[ "$passed" != 1 ]]; then
  echo "E185 failed independent seed2; skipping seeds3-5"
  touch "$E185/RESEED_DECISION_COMPLETE"
  exit 0
fi

for shard in 1 2 3; do
  found=$(find "$E185/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -ne 5 ]]; then
    bash af3_pipeline/run_worker.sh "$E185" "$shard"
  fi
  seed=$((shard + 2))
  "$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
    --experiment "$E185" \
    --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
    --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
    --shard "$shard" --out "$E185/seed${seed}_strict"
done
touch "$E185/RESEED_DECISION_COMPLETE"
echo "E185 seeds2-5 complete"
