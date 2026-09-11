#!/usr/bin/env bash
# Run the AF3-labelled 001/105 crossover screen after the 073 reversion screen.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E168="$PROJECT/experiments/e168_c4_073_reversion_ai_af3_seed1"
E164="$PROJECT/experiments/e164_c4_001_105_hit_crossover_af3_seed1"
cd "$PROJECT"

deadline=$((SECONDS + 259200))
while true; do
  complete=$(find "$E168" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  [[ "$complete" -eq 70 ]] && break
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for E168: found $complete/70 sample summaries" >&2
    exit 2
  fi
  sleep 30
done

for shard in 0 1 2 3 4 5 6 7; do
  found=$(find "$E164/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -eq 10 ]]; then
    echo "E164 shard $shard already complete"
    continue
  fi
  bash af3_pipeline/run_worker.sh "$E164" "$shard"
done

echo "E164 local seed1 screen complete"
