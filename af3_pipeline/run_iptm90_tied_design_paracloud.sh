#!/usr/bin/env bash
# Submit eight copies with unique ky-YYYYMMDD-NNN names and TASK_INDEX=0..7.
#SBATCH --partition=vip_gpu_a800_scwb671
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00

set -euo pipefail
TASK_INDEX="${TASK_INDEX:?submit with --export=ALL,TASK_INDEX=0..7}"
ROOT=/data/home/scwb671/run/wky
export OVA_TETRA_ROOT="$ROOT/ova_p3_tetramer_design"
export BINDCRAFT_MOD="$ROOT/antibody-design/pipeline_stage6_maskcap_bindcraft"
export OVA_IPTM90_EXPERIMENT="${OVA_IPTM90_EXPERIMENT:-$OVA_TETRA_ROOT/experiments/e155_c4_iptm90_crossmodel_design}"
mkdir -p "$OVA_IPTM90_EXPERIMENT/logs"
echo "$(date -u '+%F %T UTC') tied design task=$TASK_INDEX host=$(hostname)"
nvidia-smi -L
exec bash "$OVA_TETRA_ROOT/af3_pipeline/run_iptm90_tied_design_worker.sh" "$TASK_INDEX"
