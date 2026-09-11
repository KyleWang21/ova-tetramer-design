#!/usr/bin/env bash
# Priority validation of the 15-mutation E162 hit, then resume the interrupted queue.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E175="$PROJECT/experiments/e175_c4_105_rev337_356_multiseed"
cd "$PROJECT"

for shard in 0 1 2 3; do
  found=$(find "$E175/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -eq 5 ]]; then
    echo "E175 shard $shard already complete"
    continue
  fi
  bash af3_pipeline/run_worker.sh "$E175" "$shard"
done

echo "E175 priority reseed complete; resuming E162"
exec bash af3_pipeline/run_e158_then_e165_e159_e162_local.sh
