#!/usr/bin/env bash
# Fail-fast multiseed check for E174 AI-MULTIHIT-09 before E183.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E174="$PROJECT/experiments/e174_c4_073_multihit_ai_af3_seed1"
E184="$PROJECT/experiments/e184_c4_multihit09_multiseed"
PY=/root/400083/alphafold3/repo/.venv/bin/python
cd "$PROJECT"

deadline=$((SECONDS + 172800))
while true; do
  complete=$(find "$E174" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  [[ "$complete" -eq 80 ]] && break
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for E174: found $complete/80 sample summaries" >&2
    exit 2
  fi
  sleep 30
done

found=$(find "$E184/out_s0" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
if [[ "$found" -ne 5 ]]; then
  bash af3_pipeline/run_worker.sh "$E184" 0
fi
"$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
  --experiment "$E184" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
  --shard 0 --out "$E184/seed2_strict"

passed=$(awk -F '\t' '
  NR==1 {for (i=1; i<=NF; i++) if ($i=="seed1_strict_pass") column=i; next}
  NR==2 {print $column}
' "$E184/seed2_strict/seed1_strict_summary.tsv")
if [[ "$passed" != 1 ]]; then
  echo "E184 failed independent seed2; skipping seeds3-5"
  touch "$E184/RESEED_DECISION_COMPLETE"
  exit 0
fi

for shard in 1 2 3; do
  found=$(find "$E184/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -ne 5 ]]; then
    bash af3_pipeline/run_worker.sh "$E184" "$shard"
  fi
done
touch "$E184/RESEED_DECISION_COMPLETE"
echo "E184 seeds2-5 complete"
