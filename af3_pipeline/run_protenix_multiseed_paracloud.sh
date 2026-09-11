#!/usr/bin/env bash
# One Protenix-v2 shard on a Paracloud A800 SLURM compute node.
# Submit with SHARD=0..7 and a normalized --job-name=ky-YYYYMMDD-NNN.
#SBATCH --partition=vip_gpu_a800_scwb671
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00

set -uo pipefail

ROOT=/data/home/scwb671/run/wky
EXP="${OVA_PROTENIX_EXPERIMENT:-$ROOT/ova_p3_tetramer_design/experiments/e150_final20_protenix_multiseed_v1}"
MOD="$ROOT/antibody-design/pipeline_stage3_protenix_v2_refold"
SHARD="${SHARD:?submit with --export=ALL,SHARD=0..7}"
MANIFEST="$EXP/volc_manifest.tsv"
OUT="$EXP/out"

mkdir -p "$OUT" "$EXP/logs"
exec >> "$EXP/logs/protenix_multiseed_s${SHARD}.log" 2>&1
echo "$(date -u '+%F %T UTC') start shard=$SHARD host=$(hostname)"
nvidia-smi -L

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
