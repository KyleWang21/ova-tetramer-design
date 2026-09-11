#!/usr/bin/env bash
# Zhongwei one-GPU worker for the remaining v1.1 core post-screen.
# MODE is protenix or stoichiometry; OpenDDE is intentionally unsupported.
# Submit eight copies with --export=ALL,MODE=...,SHARD=0..7 and a ky-* name.
#SBATCH --partition=vip_gpu_a800_scwb671
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
set -uo pipefail

MODE="${MODE:?MODE=protenix or stoichiometry is required}"
SHARD="${SHARD:-${SLURM_ARRAY_TASK_ID:?SHARD=0..7 or a SLURM array task is required}}"
PROJECT="${ZH_PROJECT:-/data/home/scwb671/run/wky/ova_p3_tetramer_design}"
ROOT="$PROJECT/experiments/current_pending_v1_screen"
LOG="$ROOT/logs_zhongwei_${MODE}"
mkdir -p "$LOG"
exec >> "$LOG/shard_${SHARD}.log" 2>&1
echo "$(date -u '+%F %T UTC') start mode=$MODE shard=$SHARD host=$(hostname)"
nvidia-smi -L

if [[ "$MODE" == protenix ]]; then
  MANIFEST="$ROOT/zh_protenix_manifest.tsv"
  MOD="/data/home/scwb671/run/wky/antibody-design/pipeline_stage3_protenix_v2_refold"
  while IFS=$'\t' read -r row_shard candidate seed input_json; do
    row_shard="${row_shard%$'\r'}"; candidate="${candidate%$'\r'}"; seed="${seed%$'\r'}"; input_json="${input_json%$'\r'}"
    [[ "$row_shard" == shard || "$row_shard" != "$SHARD" ]] && continue
    input="$ROOT/$input_json"
    out="$ROOT/$candidate/zh_protenix_out"
    marker="$ROOT/$candidate/.ZH_PROTENIX_${seed}.DONE"
    [[ -f "$marker" ]] && continue
    mkdir -p "$out"
    echo "$(date -u '+%F %T UTC') candidate=$candidate seed=$seed"
    env CUDA_VISIBLE_DEVICES=0 LAYERNORM_TYPE=openfold USE_MSA=true \
      "$MOD/scripts/run_protenix_v2.sh" "$input" "$out" "$seed" 10 1 200
    rc=$?
    [[ "$rc" -eq 0 ]] || { echo "Protenix failed candidate=$candidate seed=$seed rc=$rc"; exit "$rc"; }
    touch "$marker"
  done < "$MANIFEST"
elif [[ "$MODE" == stoichiometry ]]; then
  MANIFEST="$ROOT/zh_stoich_manifest.tsv"
  AF3_PY="/data/home/scwb671/run/wky/alphafold3/repo/.venv/bin/python"
  AF3_RUN="/data/home/scwb671/run/wky/alphafold3/repo/run_alphafold.py"
  MODEL_DIR="/data/home/scwb671/run/wky/alphafold3/models"
  while IFS=$'\t' read -r row_shard candidate system input_json; do
    row_shard="${row_shard%$'\r'}"; candidate="${candidate%$'\r'}"; system="${system%$'\r'}"; input_json="${input_json%$'\r'}"
    [[ "$row_shard" == shard || "$row_shard" != "$SHARD" ]] && continue
    input="$ROOT/$input_json"
    out="$ROOT/$candidate/zh_stoich_out/out_s${SHARD}/$system"
    marker="$ROOT/$candidate/.ZH_STOICH_${system}.DONE"
    [[ -f "$marker" ]] && continue
    mkdir -p "$out"
    echo "$(date -u '+%F %T UTC') candidate=$candidate system=$system"
    env -u LD_LIBRARY_PATH CUDA_VISIBLE_DEVICES=0 XLA_PYTHON_CLIENT_PREALLOCATE=false \
      "$AF3_PY" "$AF3_RUN" --json_path="$input" --output_dir="$out" \
      --model_dir="$MODEL_DIR" --norun_data_pipeline --gpu_device=0
    rc=$?
    [[ "$rc" -eq 0 ]] || { echo "AF3 stoich failed candidate=$candidate system=$system rc=$rc"; exit "$rc"; }
    touch "$marker"
  done < "$MANIFEST"
else
  echo "unsupported MODE=$MODE" >&2
  exit 2
fi
echo "$(date -u '+%F %T UTC') finish mode=$MODE shard=$SHARD"
