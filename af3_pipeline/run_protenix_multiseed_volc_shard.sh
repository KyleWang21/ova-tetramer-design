#!/usr/bin/env bash
# Run one shard of the final20 Protenix-v2 multiseed screen.
set -uo pipefail

ROOT=/root/400083
EXP="${1:?usage: run_protenix_multiseed_volc_shard.sh EXPERIMENT_DIR SHARD}"
SHARD="${2:?usage: run_protenix_multiseed_volc_shard.sh EXPERIMENT_DIR SHARD}"
MOD="$ROOT/antibody-design/pipeline_stage3_protenix_v2_refold"
MANIFEST="$EXP/volc_manifest.tsv"
OUT="$EXP/out"

mkdir -p "$OUT" "$EXP/logs"
exec >> "$EXP/logs/protenix_multiseed_s${SHARD}.log" 2>&1
echo "$(date -u '+%F %T UTC') start shard=$SHARD host=$(hostname)"
nvidia-smi -L

if [ ! -e /vepfs-mlp2/c20250508/400083 ]; then
  mkdir -p /vepfs-mlp2/c20250508
  ln -s /root/400083 /vepfs-mlp2/c20250508/400083
fi

count=0
while IFS=$'\t' read -r row_shard candidate seed input_json; do
  input_json="${input_json%$'\r'}"
  [ "$row_shard" = "shard" ] && continue
  [ "$row_shard" = "$SHARD" ] || continue
  count=$((count + 1))
  [[ "$input_json" = /* ]] || input_json="$EXP/$input_json"
  echo "$(date -u '+%F %T UTC') candidate=$candidate seed=$seed"
  env CUDA_VISIBLE_DEVICES=0 LAYERNORM_TYPE=openfold USE_MSA=true \
    "$MOD/scripts/run_protenix_v2.sh" "$input_json" "$OUT" "$seed" 10 1 200
  rc=$?
  echo "$(date -u '+%F %T UTC') candidate=$candidate seed=$seed rc=$rc"
  [ "$rc" -eq 0 ] || exit "$rc"
done < "$MANIFEST"

[ "$count" -gt 0 ] || exit 2
echo "$(date -u '+%F %T UTC') finish shard=$SHARD runs=$count"
