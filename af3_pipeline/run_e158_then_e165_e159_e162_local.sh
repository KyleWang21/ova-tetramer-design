#!/usr/bin/env bash
# Preserve one-GPU execution while prioritizing the stronger geometry-clean
# E158 hit (073), then complete the earlier 105 reseed and two seed-1 libraries.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E158="$PROJECT/experiments/e158_c4_iptm90_joint_surrogate_af3_seed1"
E165="$PROJECT/experiments/e165_c4_iptm90_073_multiseed"
E159="$PROJECT/experiments/e159_c4_iptm90_105_multiseed"
E162="$PROJECT/experiments/e162_c4_105_reversion_ai_af3_seed1"
cd "$PROJECT"

while pgrep -f '[r]un_alphafold.py.*e158_c4_iptm90_joint_surrogate_af3_seed1' >/dev/null; do
  sleep 30
done

complete=$(find "$E158/out_s0" -path '*/seed-*/*_summary_confidences.json' | wc -l)
if [[ "$complete" -ne 80 ]]; then
  echo "E158 shard-0 incomplete: expected 80 sample summaries, found $complete" >&2
  exit 2
fi

for shard in 0 1 2 3; do
  found=$(find "$E165/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -eq 5 ]]; then
    echo "E165 seed $((shard + 2)) already complete"
    continue
  fi
  bash af3_pipeline/run_worker.sh "$E165" "$shard"
done

for shard in 0 1 2 3; do
  found=$(find "$E159/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -eq 5 ]]; then
    echo "E159 seed $((shard + 2)) already complete"
    continue
  fi
  bash af3_pipeline/run_worker.sh "$E159" "$shard"
done

for shard in 0 1 2 3 4 5 6 7; do
  found=$(find "$E162/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -eq 10 ]]; then
    echo "E162 shard $shard already complete"
    continue
  fi
  bash af3_pipeline/run_worker.sh "$E162" "$shard"
done

echo "E165, E159 and E162 local queue complete"
