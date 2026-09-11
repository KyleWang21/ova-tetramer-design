#!/usr/bin/env bash
# Run the formal tied-ProteinMPNN sequence library after all queued E174 hits.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E189="$PROJECT/experiments/e189_c4_multihit15_multiseed"
E177="$PROJECT/experiments/e177_c4_tied_proteinmpnn_af3_seed1"
VOLC_MARKER="$E177/VOLC_AF3_SHARDS_SUBMITTED"
cd "$PROJECT"

deadline=$((SECONDS + 259200))
while true; do
  [[ -f "$E189/RESEED_DECISION_COMPLETE" ]] && break
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for E189 reseed decision" >&2
    exit 2
  fi
  sleep 30
done

if [[ -f "$VOLC_MARKER" ]]; then
  while true; do
    complete=$(find "$E177" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
    [[ "$complete" -eq 160 ]] && exit 0
    if (( SECONDS >= deadline )); then
      echo "Volc AF3 wait timed out at $complete/160; falling back to local shards" >&2
      break
    fi
    sleep 30
  done
fi

for shard in 0 1 2 3 4 5 6 7; do
  found=$(find "$E177/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -eq 20 ]]; then
    echo "E177 shard $shard already complete"
    continue
  fi
  bash af3_pipeline/run_worker.sh "$E177" "$shard"
done

echo "E177 tied-ProteinMPNN AF3 seed-1 screen complete"
