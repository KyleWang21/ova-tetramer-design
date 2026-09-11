#!/usr/bin/env bash
# Run five independent Protenix-v2 seeds for 073 after the queued local AF3 work.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E171="$PROJECT/experiments/e171_c4_073_stoichiometry_v1"
E170="$PROJECT/experiments/e170_c4_073_protenix_v2_multiseed"
MOD=/root/400083/antibody-design/pipeline_stage3_protenix_v2_refold
INPUT="$E170/base/input/ova_iptm80_ova_body_c4_073_ms_protenix.json"
OUT="$E170/out_local"
cd "$PROJECT"

deadline=$((SECONDS + 259200))
while true; do
  complete=$(find "$E171" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  [[ "$complete" -eq 75 ]] && break
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for E171: found $complete/75 sample summaries" >&2
    exit 2
  fi
  sleep 30
done

mkdir -p "$OUT" "$E170/logs"
for seed in 11 21 31 41 51; do
  found=$(find "$OUT" -path "*/seed_${seed}/*_summary_confidence_sample_0.json" 2>/dev/null | wc -l || true)
  if [[ "$found" -eq 1 ]]; then
    echo "E170 seed $seed already complete"
    continue
  fi
  env CUDA_VISIBLE_DEVICES=0 LAYERNORM_TYPE=openfold USE_MSA=true \
    "$MOD/scripts/run_protenix_v2.sh" "$INPUT" "$OUT" "$seed" 10 1 200
done

echo "E170 073 five-seed Protenix-v2 screen complete"
