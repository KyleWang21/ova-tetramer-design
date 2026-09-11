#!/usr/bin/env bash
# Run direct AF3-medoid tied-ProteinMPNN proposals after higher-priority screens.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E164="$PROJECT/experiments/e164_c4_001_105_hit_crossover_af3_seed1"
E181="$PROJECT/experiments/e181_c4_073_af3medoid_mpnn_af3_seed1"
cd "$PROJECT"

deadline=$((SECONDS + 345600))
while true; do
  complete=$(find "$E164" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  [[ "$complete" -eq 80 ]] && break
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for E164: found $complete/80 sample summaries" >&2
    exit 2
  fi
  sleep 30
done

for shard in 0 1 2 3 4 5 6 7; do
  found=$(find "$E181/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -eq 10 ]]; then
    echo "E181 shard $shard already complete"
    continue
  fi
  bash af3_pipeline/run_worker.sh "$E181" "$shard"
done

echo "E181 AF3-medoid tied-ProteinMPNN seed-1 screen complete"
