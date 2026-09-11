#!/usr/bin/env bash
# Run one shard of the final OVA Protenix-v2 screen on a Volcengine A100.
set -uo pipefail

ROOT=/root/400083
PROJECT="$ROOT/ova_p3_tetramer_design"
EXP="${1:?usage: run_protenix_volc_shard.sh EXPERIMENT_DIR SHARD}"
SHARD="${2:?usage: run_protenix_volc_shard.sh EXPERIMENT_DIR SHARD}"
MOD="$ROOT/antibody-design/pipeline_stage3_protenix_v2_refold"
MANIFEST="$EXP/volc_manifest.tsv"
OUT="$EXP/out_volc"

mkdir -p "$OUT" "$EXP/logs"
exec >> "$EXP/logs/protenix_volc_s${SHARD}.log" 2>&1
echo "$(date -u '+%F %T UTC') start shard=$SHARD host=$(hostname)"
nvidia-smi -L

if [ ! -e /vepfs-mlp2/c20250508/400083 ]; then
  mkdir -p /vepfs-mlp2/c20250508
  ln -s /root/400083 /vepfs-mlp2/c20250508/400083
fi

count=0
while IFS=$'\t' read -r row_shard candidate input_json; do
  [ "$row_shard" = "shard" ] && continue
  [ "$row_shard" = "$SHARD" ] || continue
  count=$((count + 1))
  echo "$(date -u '+%F %T UTC') candidate=$candidate input=$input_json"
  env CUDA_VISIBLE_DEVICES=0 LAYERNORM_TYPE=openfold USE_MSA=true \
    "$MOD/scripts/run_protenix_v2.sh" "$input_json" "$OUT" 31 10 1 200
  rc=$?
  echo "$(date -u '+%F %T UTC') candidate=$candidate rc=$rc"
  [ "$rc" -eq 0 ] || exit "$rc"
done < "$MANIFEST"

if [ "$count" -eq 0 ]; then
  echo "no candidates for shard=$SHARD"
  exit 2
fi
echo "$(date -u '+%F %T UTC') finish shard=$SHARD candidates=$count rc=0"
