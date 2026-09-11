#!/usr/bin/env bash
# Queue the four independent reseeds of the first geometry-clean E158 hit on
# the same local A100, but only after E158 shard-0 completed all 16 systems.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E158="$PROJECT/experiments/e158_c4_iptm90_joint_surrogate_af3_seed1"
E159="$PROJECT/experiments/e159_c4_iptm90_105_multiseed"
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
  expected=$((shard + 2))
  found=$(find "$E159/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$found" -eq 5 ]]; then
    echo "E159 seed $expected already complete"
    continue
  fi
  bash af3_pipeline/run_worker.sh "$E159" "$shard"
done

echo "E159 seeds2-5 local queue complete"
