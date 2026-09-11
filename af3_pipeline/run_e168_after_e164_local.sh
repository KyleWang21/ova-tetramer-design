#!/usr/bin/env bash
# Run the AF3-interface-derived 073 reversion library after E178.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E178="$PROJECT/experiments/e178_c4_073_af3ensemble_ai_af3_seed1"
E168="$PROJECT/experiments/e168_c4_073_reversion_ai_af3_seed1"
cd "$PROJECT"

deadline=$((SECONDS + 86400))
while true; do
  complete=$(find "$E178" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  [[ "$complete" -eq 80 ]] && break
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for E178: found $complete/80 sample summaries" >&2
    exit 2
  fi
  sleep 30
done

for shard in 0 1 2 3 4 5 6; do
  found=$(find "$E168/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -eq 10 ]]; then
    echo "E168 shard $shard already complete"
    continue
  fi
  bash af3_pipeline/run_worker.sh "$E168" "$shard"
done

echo "E168 local seed1 screen complete"
