#!/usr/bin/env bash
# One official OpenDDE-v1 shard on a Paracloud A800. Submit with a normalized
# --job-name=ky-YYYYMMDD-NNN and SHARD=0..7.
#SBATCH --partition=vip_gpu_a800_scwb671
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00

set -euo pipefail

ROOT=/data/home/scwb671/run/wky
PROJECT="$ROOT/ova_p3_tetramer_design"
EXP="${EXP:-$PROJECT/experiments/e152_final20_opendde_v1}"
SHARD="${SHARD:?submit with --export=ALL,SHARD=0..7}"
export OPENDDE_ROOT_DIR="$ROOT/OpenDDE_data"
export OPENDDE_BIN="$ROOT/antibody-design/pipeline_stage3_protenix_v2_refold/venv/bin/opendde"
# OpenDDE and Protenix both ship a top-level Python package named `runner`.
# This environment is reused from Protenix, so force the released OpenDDE
# source tree to the front or the CLI silently imports Protenix's runner.
export PYTHONPATH="$ROOT/OpenDDE${PYTHONPATH:+:$PYTHONPATH}"
export MANIFEST="${MANIFEST:-$EXP/manifest.tsv}"

mkdir -p "$EXP/logs"
echo "$(date -u '+%F %T UTC') OpenDDE shard=$SHARD host=$(hostname) manifest=$MANIFEST"
nvidia-smi -L
exec bash "$PROJECT/af3_pipeline/run_opendde_c4_shard.sh" "$EXP" "$SHARD"
