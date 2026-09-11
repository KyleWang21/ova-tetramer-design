#!/usr/bin/env bash
# Run the 073 C1/C2/C3/C5/C6 controls after all queued AF3 design screens.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E181="$PROJECT/experiments/e181_c4_073_af3medoid_mpnn_af3_seed1"
E171="$PROJECT/experiments/e171_c4_073_stoichiometry_v1"
cd "$PROJECT"

deadline=$((SECONDS + 172800))
while true; do
  complete=$(find "$E181" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  [[ "$complete" -eq 80 ]] && break
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for E181: found $complete/80 sample summaries" >&2
    exit 2
  fi
  sleep 30
done

for shard in 0 1 2 3 4; do
  expected=15
  found=$(find "$E171/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -eq "$expected" ]]; then
    echo "E171 shard $shard already complete"
    continue
  fi
  bash af3_pipeline/run_worker.sh "$E171" "$shard"
done

echo "E171 073 stoichiometry screen complete"
