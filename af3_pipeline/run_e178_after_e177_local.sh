#!/usr/bin/env bash
# Run AF3-coordinate tied-AF2 proposals after the ProteinMPNN screen.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E177="$PROJECT/experiments/e177_c4_tied_proteinmpnn_af3_seed1"
E178="$PROJECT/experiments/e178_c4_073_af3ensemble_ai_af3_seed1"
FROZEN="$E178/CANDIDATES_FROZEN"
VOLC_MARKER="$E178/VOLC_AF3_SHARDS_SUBMITTED"
cd "$PROJECT"

deadline=$((SECONDS + 259200))
while true; do
  complete=$(find "$E177" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  [[ "$complete" -eq 160 ]] && break
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for E177: found $complete/160 sample summaries" >&2
    exit 2
  fi
  sleep 30
done

while [[ ! -f "$FROZEN" ]]; do
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for refreshed E178 candidate freeze marker" >&2
    exit 2
  fi
  sleep 30
done

if [[ -f "$VOLC_MARKER" ]]; then
  while true; do
    complete=$(find "$E178" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
    [[ "$complete" -eq 80 ]] && exit 0
    if (( SECONDS >= deadline )); then
      echo "Volc E178 wait timed out at $complete/80; falling back to local shards" >&2
      break
    fi
    sleep 30
  done
fi

for shard in 0 1 2 3 4 5 6 7; do
  found=$(find "$E178/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -eq 10 ]]; then
    echo "E178 shard $shard already complete"
    continue
  fi
  bash af3_pipeline/run_worker.sh "$E178" "$shard"
done

echo "E178 AF3-coordinate tied-AF2 seed-1 screen complete"
