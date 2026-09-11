#!/usr/bin/env bash
# Run the AF3-ranked interface-reversion library after all four E159 reseeds.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E159="$PROJECT/experiments/e159_c4_iptm90_105_multiseed"
E162="$PROJECT/experiments/e162_c4_105_reversion_ai_af3_seed1"
cd "$PROJECT"

deadline=$((SECONDS + 43200))
while true; do
  complete=$(find "$E159" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  [[ "$complete" -eq 20 ]] && break
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for E159: found $complete/20 sample summaries" >&2
    exit 2
  fi
  sleep 30
done

for shard in 0 1 2 3 4 5 6 7; do
  found=$(find "$E162/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -eq 10 ]]; then
    echo "E162 shard $shard already complete"
    continue
  fi
  bash af3_pipeline/run_worker.sh "$E162" "$shard"
done

echo "E162 local seed1 screen complete"
