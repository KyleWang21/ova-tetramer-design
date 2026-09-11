#!/usr/bin/env bash
# Run the six-parent AF3-labelled allele library after E180's reseed decision.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E180="$PROJECT/experiments/e180_c4_105_rev155_337_356_multiseed"
E174="$PROJECT/experiments/e174_c4_073_multihit_ai_af3_seed1"
cd "$PROJECT"

deadline=$((SECONDS + 172800))
while true; do
  [[ -f "$E180/RESEED_DECISION_COMPLETE" ]] && break
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for E180 reseed decision" >&2
    exit 2
  fi
  sleep 30
done

for shard in 0 1 2 3 4 5 6 7; do
  found=$(find "$E174/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -eq 10 ]]; then
    echo "E174 shard $shard already complete"
    continue
  fi
  bash af3_pipeline/run_worker.sh "$E174" "$shard"
done

echo "E174 multi-hit AI seed1 screen complete"
